#include "qRoadMarking.h"

#include <ccHObjectCaster.h>
#include <ccMainAppInterface.h>
#include <ccPointCloud.h>

#include <ScalarField.h>

#include <QAction>
#include <QMessageBox>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <queue>

qRoadMarking::qRoadMarking(QObject* parent)
    : QObject(parent)
    , ccStdPluginInterface("../info.json")
{
    m_action = new QAction(tr("Extract Road Markings (Pan 2019)"), this);
    m_action->setToolTip(tr("Pan-style road-marking extraction (OTSU intensity + component filtering)"));

    connect(m_action, &QAction::triggered, this, &qRoadMarking::doAction);
}

void qRoadMarking::onNewSelection(const ccHObject::Container& selectedEntities)
{
    m_hasSelection = false;
    for (ccHObject* entity : selectedEntities)
    {
        if (ccHObjectCaster::ToPointCloud(entity))
        {
            m_hasSelection = true;
            break;
        }
    }

    if (m_action)
    {
        m_action->setEnabled(m_hasSelection);
    }
}

QList<QAction*> qRoadMarking::getActions()
{
    if (m_action)
    {
        m_action->setEnabled(m_hasSelection);
        return {m_action};
    }

    return {};
}

ccPointCloud* qRoadMarking::getSelectedPointCloud() const
{
    if (!m_app)
    {
        return nullptr;
    }

    const ccHObject::Container& selected = m_app->getSelectedEntities();
    for (ccHObject* entity : selected)
    {
        if (ccPointCloud* cloud = ccHObjectCaster::ToPointCloud(entity))
        {
            return cloud;
        }
    }

    return nullptr;
}

void qRoadMarking::doAction()
{
    ccPointCloud* inputCloud = getSelectedPointCloud();
    if (!inputCloud)
    {
        QMessageBox::warning(nullptr, tr("Road marking"), tr("Please select a point cloud first."));
        return;
    }

    const int sfIndex = inputCloud->getCurrentDisplayedScalarFieldIndex();
    if (sfIndex < 0)
    {
        QMessageBox::warning(nullptr,
                             tr("Road marking"),
                             tr("No scalar field selected. Please select intensity as displayed scalar field."));
        return;
    }

    const CCCoreLib::ScalarField* intensitySf = inputCloud->getScalarField(sfIndex);
    if (!intensitySf)
    {
        QMessageBox::warning(nullptr, tr("Road marking"), tr("Unable to access current scalar field."));
        return;
    }

    PanParams params;
    std::vector<unsigned char> mask = PanRoadMarkingAdapter::extractMask(inputCloud, intensitySf, params);
    PanRoadMarkingAdapter::sorFilter(inputCloud, mask, params.meanK, params.sorStdMul);

    auto* outputCloud = new ccPointCloud(QString("%1_road_markings").arg(inputCloud->getName()));
    if (!outputCloud->reserveThePointsTable(inputCloud->size()))
    {
        QMessageBox::critical(nullptr, tr("Road marking"), tr("Not enough memory for output cloud."));
        delete outputCloud;
        return;
    }

    for (unsigned i = 0; i < inputCloud->size(); ++i)
    {
        if (i < mask.size() && mask[i] != 0)
        {
            outputCloud->addPoint(*inputCloud->getPoint(i));
        }
    }

    outputCloud->shrinkToFit();

    m_app->addToDB(outputCloud);
    m_app->dispToConsole(tr("RoadMarking: %1 / %2 points extracted")
                             .arg(static_cast<qulonglong>(outputCloud->size()))
                             .arg(static_cast<qulonglong>(inputCloud->size())),
                         ccMainAppInterface::STD_CONSOLE_MESSAGE);
    m_app->refreshAll();
}

std::vector<unsigned char> PanRoadMarkingAdapter::extractMask(const ccPointCloud* cloud,
                                                              const CCCoreLib::ScalarField* intensitySf,
                                                              const PanParams& params)
{
    if (!cloud || !intensitySf)
    {
        return {};
    }

    const unsigned count = cloud->size();
    if (count == 0)
    {
        return {};
    }

    std::vector<unsigned char> mask(count, 0);

    ScalarType maxIntensity = -FLT_MAX;
    for (unsigned i = 0; i < count; ++i)
    {
        const ScalarType v = intensitySf->getValue(i);
        if (CCCoreLib::ScalarField::ValidValue(v))
        {
            maxIntensity = std::max(maxIntensity, v);
        }
    }

    if (maxIntensity <= 0)
    {
        return mask;
    }

    const int bins = std::max(2, params.histogramBins);
    std::vector<int> histogram(static_cast<size_t>(bins), 0);

    for (unsigned i = 0; i < count; ++i)
    {
        const ScalarType v = intensitySf->getValue(i);
        if (!CCCoreLib::ScalarField::ValidValue(v) || v < 0)
        {
            continue;
        }

        const int bin = std::clamp(static_cast<int>((bins - 1) * v / maxIntensity), 0, bins - 1);
        histogram[static_cast<size_t>(bin)] += 1;
    }

    int thresholdBin = bins - 1;
    if (params.useOtsuThreshold)
    {
        thresholdBin = computeOtsuThresholdBin(histogram, static_cast<int>(count));
    }
    else
    {
        thresholdBin = std::clamp(static_cast<int>((bins - 1) * params.manualIntensityThreshold / maxIntensity), 0, bins - 1);
    }

    for (unsigned i = 0; i < count; ++i)
    {
        const ScalarType v = intensitySf->getValue(i);
        if (!CCCoreLib::ScalarField::ValidValue(v) || v < 0)
        {
            continue;
        }

        const int bin = std::clamp(static_cast<int>((bins - 1) * v / maxIntensity), 0, bins - 1);
        if (bin > thresholdBin)
        {
            mask[i] = 1;
        }
    }

    return mask;
}

int PanRoadMarkingAdapter::computeOtsuThresholdBin(const std::vector<int>& histogram, int pointCount)
{
    if (histogram.empty() || pointCount <= 0)
    {
        return 0;
    }

    const int bins = static_cast<int>(histogram.size());
    std::vector<double> prob(static_cast<size_t>(bins), 0.0);
    for (int i = 0; i < bins; ++i)
    {
        prob[static_cast<size_t>(i)] = static_cast<double>(histogram[static_cast<size_t>(i)]) / pointCount;
    }

    double best = -1.0;
    int threshold = 0;

    for (int k = 0; k < bins; ++k)
    {
        double w0 = 0.0;
        double w1 = 0.0;
        double u0acc = 0.0;
        double u1acc = 0.0;

        for (int t = 0; t < bins; ++t)
        {
            if (t <= k)
            {
                w0 += prob[static_cast<size_t>(t)];
                u0acc += t * prob[static_cast<size_t>(t)];
            }
            else
            {
                w1 += prob[static_cast<size_t>(t)];
                u1acc += t * prob[static_cast<size_t>(t)];
            }
        }

        if (w0 <= 0.0 || w1 <= 0.0)
        {
            continue;
        }

        const double u0 = u0acc / w0;
        const double u1 = u1acc / w1;
        const double u = u0acc + u1acc;
        const double delta = w0 * (u0 - u) * (u0 - u) + w1 * (u1 - u) * (u1 - u);

        if (delta > best)
        {
            best = delta;
            threshold = k;
        }
    }

    return threshold;
}

float PanRoadMarkingAdapter::sqDist(const ccPointCloud* cloud, unsigned i, unsigned j)
{
    const CCVector3* pi = cloud->getPoint(i);
    const CCVector3* pj = cloud->getPoint(j);
    const float dx = pi->x - pj->x;
    const float dy = pi->y - pj->y;
    const float dz = pi->z - pj->z;
    return dx * dx + dy * dy + dz * dz;
}

void PanRoadMarkingAdapter::sorFilter(const ccPointCloud* cloud,
                                      std::vector<unsigned char>& mask,
                                      int meanK,
                                      float stdMul)
{
    if (!cloud || mask.empty() || meanK <= 1)
    {
        return;
    }

    std::vector<unsigned> indices;
    indices.reserve(mask.size());
    for (unsigned i = 0; i < mask.size(); ++i)
    {
        if (mask[i] != 0)
        {
            indices.push_back(i);
        }
    }

    if (indices.size() <= static_cast<size_t>(meanK))
    {
        return;
    }

    std::vector<float> meanD(indices.size(), 0.0f);
    std::vector<float> knn;
    knn.reserve(indices.size());

    for (size_t ii = 0; ii < indices.size(); ++ii)
    {
        knn.clear();
        const unsigned i = indices[ii];

        for (size_t jj = 0; jj < indices.size(); ++jj)
        {
            if (ii == jj)
            {
                continue;
            }
            knn.push_back(sqDist(cloud, i, indices[jj]));
        }

        const size_t k = std::min(static_cast<size_t>(meanK), knn.size());
        std::nth_element(knn.begin(), knn.begin() + (k - 1), knn.end());

        float sum = 0.0f;
        for (size_t t = 0; t < k; ++t)
        {
            sum += std::sqrt(knn[t]);
        }
        meanD[ii] = sum / static_cast<float>(k);
    }

    float mu = 0.0f;
    for (float d : meanD)
    {
        mu += d;
    }
    mu /= static_cast<float>(meanD.size());

    float var = 0.0f;
    for (float d : meanD)
    {
        const float e = d - mu;
        var += e * e;
    }
    var /= static_cast<float>(meanD.size());
    const float sigma = std::sqrt(var);
    const float maxDist = mu + stdMul * sigma;

    for (size_t ii = 0; ii < indices.size(); ++ii)
    {
        if (meanD[ii] > maxDist)
        {
            mask[indices[ii]] = 0;
        }
    }
}
