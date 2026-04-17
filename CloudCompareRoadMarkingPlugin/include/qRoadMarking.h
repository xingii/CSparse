#pragma once

#include <ccStdPluginInterface.h>

#include <QAction>
#include <QObject>
#include <vector>

class ccPointCloud;

namespace CCCoreLib
{
class ScalarField;
}

struct PanParams
{
    // 与 Pan et al. 2019 代码中的强度直方图思路一致：优先使用 OTSU 自动阈值
    bool useOtsuThreshold = true;
    float manualIntensityThreshold = 35.0f;
    int histogramBins = 256;

    // Pan cloudFilter 参数风格：OTSU + SOR
    int meanK = 12;
    float sorStdMul = 1.0f;
};

class PanRoadMarkingAdapter
{
public:
    // 基于 Pan 代码中 cloudFilter 的核心思路：
    // 1) intensity histogram + OTSU
    // 2) threshold segmentation
    static std::vector<unsigned char> extractMask(const ccPointCloud* cloud,
                                                  const CCCoreLib::ScalarField* intensitySf,
                                                  const PanParams& params);

    static void sorFilter(const ccPointCloud* cloud,
                          std::vector<unsigned char>& mask,
                          int meanK,
                          float stdMul);

private:
    static int computeOtsuThresholdBin(const std::vector<int>& histogram, int pointCount);
    static float sqDist(const ccPointCloud* cloud, unsigned i, unsigned j);
};

class qRoadMarking final : public QObject, public ccStdPluginInterface
{
    Q_OBJECT
    Q_INTERFACES(ccPluginInterface)
    Q_PLUGIN_METADATA(IID "cccorp.cloudcompare.plugin.qRoadMarking" FILE "../info.json")

public:
    explicit qRoadMarking(QObject* parent = nullptr);
    ~qRoadMarking() override = default;

    // ccPluginInterface
    void onNewSelection(const ccHObject::Container& selectedEntities) override;

    // ccStdPluginInterface
    QList<QAction*> getActions() override;

private slots:
    void doAction();

private:
    ccPointCloud* getSelectedPointCloud() const;

private:
    QAction* m_action = nullptr;
    bool m_hasSelection = false;
};

