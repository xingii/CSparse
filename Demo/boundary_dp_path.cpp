#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

// --------------------------
// 输入数据结构（与题述保持一致）
// --------------------------
struct CCVector2
{
    float x = 0.0f;
    float y = 0.0f;
};

struct CCVector3
{
    float x = 0.0f;
    float y = 0.0f;
    float z = 0.0f;
};

typedef struct tagSLPoint
{
    int sid;      // 数据点所属扫描线 id
    float offset; // 点到轨迹线距离
    CCVector3 pt; // 点坐标
} SLPoint;

namespace boundary_dp
{

struct Params
{
    int max_gap = 10;                 // 线段连接最大扫描线间隔
    float angle_threshold_deg = 45.f; // 条件(1)中的阈值 a
    float distance_threshold = 0.50f; // 条件(2)(3)中的阈值 d

    float w_data = 1.0f;   // 数据代价权重
    float w_smooth = 1.0f; // 平滑代价权重

    float w_length = 1.0f;   // 数据代价中长度项权重（越长越优）
    float w_parallel = 1.0f; // 数据代价中平行项权重（越平行越优）
};

struct SegmentGeom
{
    int id = -1;
    std::vector<int> point_indices;

    int start_sid = 0;
    int end_sid = 0;

    CCVector2 start_xy{};
    CCVector2 end_xy{};
    CCVector2 dir{};      // 单位方向向量
    float length = 0.0f;  // 2D 长度
    float data_cost = 0.f;
};

struct Edge
{
    int from = -1;
    int to = -1;
    float cost = 0.f; // 仅平滑代价，节点代价单独加
};

struct PathResult
{
    std::vector<int> segment_path; // 最优路径上的线段 id 序列
    float total_cost = std::numeric_limits<float>::infinity();
    bool success = false;
};

static inline CCVector2 to2D(const CCVector3 &p)
{
    return CCVector2{p.x, p.y};
}

static inline CCVector2 sub(const CCVector2 &a, const CCVector2 &b)
{
    return CCVector2{a.x - b.x, a.y - b.y};
}

static inline float dot(const CCVector2 &a, const CCVector2 &b)
{
    return a.x * b.x + a.y * b.y;
}

static inline float cross(const CCVector2 &a, const CCVector2 &b)
{
    return a.x * b.y - a.y * b.x;
}

static inline float norm(const CCVector2 &v)
{
    return std::sqrt(dot(v, v));
}

static inline CCVector2 normalize(const CCVector2 &v)
{
    float n = norm(v);
    if (n < 1e-6f)
    {
        return CCVector2{1.f, 0.f};
    }
    return CCVector2{v.x / n, v.y / n};
}

static inline float angleDeg(const CCVector2 &u, const CCVector2 &v)
{
    float c = dot(normalize(u), normalize(v));
    c = std::max(-1.0f, std::min(1.0f, c));
    return std::acos(c) * 180.0f / static_cast<float>(M_PI);
}

// 计算点 P 到直线 AB 的垂直距离（非线段距离）
static float pointLineDistance(const CCVector2 &p, const CCVector2 &a, const CCVector2 &b)
{
    CCVector2 ab = sub(b, a);
    float abn = norm(ab);
    if (abn < 1e-6f)
    {
        return norm(sub(p, a));
    }
    return std::fabs(cross(sub(p, a), ab)) / abn;
}

// 点 P 在有向直线 AB 的哪一侧：>0 左侧，<0 右侧，=0 共线
static float sideOfLine(const CCVector2 &p, const CCVector2 &a, const CCVector2 &b)
{
    return cross(sub(b, a), sub(p, a));
}

static CCVector2 queryTrajectoryTangent(const std::map<int, CCVector2> &trajectory_pts, int sid)
{
    if (trajectory_pts.empty())
    {
        return CCVector2{1.f, 0.f};
    }

    auto it = trajectory_pts.lower_bound(sid);
    if (it == trajectory_pts.begin())
    {
        auto it2 = std::next(it);
        if (it2 == trajectory_pts.end())
        {
            return CCVector2{1.f, 0.f};
        }
        return normalize(sub(it2->second, it->second));
    }
    if (it == trajectory_pts.end())
    {
        auto last = std::prev(it);
        if (last == trajectory_pts.begin())
        {
            return CCVector2{1.f, 0.f};
        }
        auto prev = std::prev(last);
        return normalize(sub(last->second, prev->second));
    }

    auto prev = std::prev(it);
    return normalize(sub(it->second, prev->second));
}

static SegmentGeom buildSegmentGeom(
    int seg_id,
    const std::vector<int> &indices,
    const std::vector<SLPoint> &pts,
    const std::map<int, CCVector2> &trajectory_pts,
    const Params &params)
{
    if (indices.empty())
    {
        throw std::runtime_error("segment contains no point indices");
    }

    SegmentGeom g;
    g.id = seg_id;
    g.point_indices = indices;

    // 按 sid 排序，首尾点用于构造线段
    std::vector<int> sorted = indices;
    std::sort(sorted.begin(), sorted.end(), [&](int a, int b) {
        return pts[a].sid < pts[b].sid;
    });

    int i0 = sorted.front();
    int i1 = sorted.back();
    g.start_sid = pts[i0].sid;
    g.end_sid = pts[i1].sid;

    g.start_xy = to2D(pts[i0].pt);
    g.end_xy = to2D(pts[i1].pt);

    CCVector2 v = sub(g.end_xy, g.start_xy);
    g.length = norm(v);
    g.dir = normalize(v);

    // 数据代价：长度越长越好（代价更小），与轨迹越平行越好（代价更小）
    CCVector2 traj_dir = queryTrajectoryTangent(trajectory_pts, (g.start_sid + g.end_sid) / 2);
    float parallel = std::fabs(dot(g.dir, traj_dir)); // [0,1]

    float length_term = params.w_length / std::max(g.length, 1e-3f);
    float parallel_term = params.w_parallel * (1.0f - parallel);

    g.data_cost = length_term + parallel_term;
    return g;
}

// 拓扑平滑检查 + 平滑代价
static bool smoothCostByTopology(
    const SegmentGeom &A,
    const SegmentGeom &B,
    const Params &params,
    float &smooth_cost)
{
    smooth_cost = 0.f;

    // 要求按扫描线方向前进：B 必须在 A 后面
    if (B.start_sid <= A.end_sid)
    {
        return false;
    }

    // 先检查转角阈值 (1)
    float angle = angleDeg(A.dir, B.dir);
    if (angle > params.angle_threshold_deg)
    {
        return false;
    }

    // 按题意选择较长 L1 和较短 L2
    const SegmentGeom *L1 = &A;
    const SegmentGeom *L2 = &B;
    if (A.length < B.length)
    {
        L1 = &B;
        L2 = &A;
    }

    float d_head = pointLineDistance(L2->start_xy, L1->start_xy, L1->end_xy);

    // (2) L2 首点到 L1 的垂直距离不大于 d，则可连接
    if (d_head <= params.distance_threshold)
    {
        smooth_cost = angle / std::max(params.angle_threshold_deg, 1e-3f);
        return true;
    }

    // (3) d_head > d，比较“同侧性”
    float side_head = sideOfLine(L2->start_xy, L1->start_xy, L1->end_xy);
    float side_end = sideOfLine(L2->end_xy, L1->end_xy, L2->start_xy);

    bool same_side = (side_head == 0.0f || side_end == 0.0f) || (side_head * side_end > 0.0f);
    if (!same_side)
    {
        return false;
    }

    smooth_cost = angle / std::max(params.angle_threshold_deg, 1e-3f);
    return true;
}

static bool overlapOnTrajectory(const SegmentGeom &A, const SegmentGeom &B)
{
    // 扫描线区间重叠则禁止连接
    return !(A.end_sid < B.start_sid || B.end_sid < A.start_sid);
}

static std::vector<Edge> buildEdges(const std::vector<SegmentGeom> &segs, const Params &params)
{
    std::vector<Edge> edges;
    const int n = static_cast<int>(segs.size());

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (i == j)
            {
                continue;
            }

            const SegmentGeom &A = segs[i];
            const SegmentGeom &B = segs[j];

            // 保证前后顺序，且间隔不超过 max_gap
            if (B.start_sid <= A.end_sid)
            {
                continue;
            }
            if (B.start_sid - A.end_sid > params.max_gap)
            {
                continue;
            }

            // 轨迹区间不能重叠
            if (overlapOnTrajectory(A, B))
            {
                continue;
            }

            float smooth = 0.f;
            if (!smoothCostByTopology(A, B, params, smooth))
            {
                continue;
            }

            edges.push_back(Edge{i, j, params.w_smooth * smooth});
        }
    }

    return edges;
}

PathResult findBestBoundaryPath(
    const std::map<int, CCVector2> &trajectory_pts,
    const std::vector<SLPoint> &pts,
    const std::vector<std::vector<int>> &segments,
    const Params &params)
{
    PathResult result;
    if (segments.empty())
    {
        return result;
    }

    std::vector<SegmentGeom> geoms;
    geoms.reserve(segments.size());
    for (int i = 0; i < static_cast<int>(segments.size()); ++i)
    {
        geoms.push_back(buildSegmentGeom(i, segments[i], pts, trajectory_pts, params));
    }

    std::vector<Edge> edges = buildEdges(geoms, params);

    // 为 DP 构建邻接表（图是 DAG：按 sid 单调递增）
    const int n = static_cast<int>(geoms.size());
    std::vector<std::vector<std::pair<int, float>>> adj(n);
    std::vector<int> indeg(n, 0);
    for (const Edge &e : edges)
    {
        adj[e.from].push_back({e.to, e.cost});
        indeg[e.to]++;
    }

    // 拓扑排序（Kahn）
    std::vector<int> topo;
    topo.reserve(n);
    std::vector<int> q;
    q.reserve(n);
    for (int i = 0; i < n; ++i)
    {
        if (indeg[i] == 0)
        {
            q.push_back(i);
        }
    }
    for (size_t h = 0; h < q.size(); ++h)
    {
        int u = q[h];
        topo.push_back(u);
        for (const auto &vw : adj[u])
        {
            int v = vw.first;
            if (--indeg[v] == 0)
            {
                q.push_back(v);
            }
        }
    }

    // s->i 的初值：选择任一线段作为起点，代价为数据代价
    std::vector<float> dp(n, std::numeric_limits<float>::infinity());
    std::vector<int> pre(n, -1);
    for (int i = 0; i < n; ++i)
    {
        dp[i] = params.w_data * geoms[i].data_cost;
    }

    // DAG 上松弛：dp[v] = min(dp[u] + smooth(u,v) + data(v))
    for (int u : topo)
    {
        if (!std::isfinite(dp[u]))
        {
            continue;
        }
        for (const auto &vw : adj[u])
        {
            int v = vw.first;
            float trans = vw.second;
            float cand = dp[u] + trans + params.w_data * geoms[v].data_cost;
            if (cand < dp[v])
            {
                dp[v] = cand;
                pre[v] = u;
            }
        }
    }

    // i->t：这里设 0 代价，选择最小 dp[i] 作为终点
    int best = -1;
    float best_cost = std::numeric_limits<float>::infinity();
    for (int i = 0; i < n; ++i)
    {
        if (dp[i] < best_cost)
        {
            best_cost = dp[i];
            best = i;
        }
    }

    if (best < 0 || !std::isfinite(best_cost))
    {
        return result;
    }

    std::vector<int> rev;
    for (int cur = best; cur >= 0; cur = pre[cur])
    {
        rev.push_back(geoms[cur].id);
    }
    std::reverse(rev.begin(), rev.end());

    result.segment_path = std::move(rev);
    result.total_cost = best_cost;
    result.success = true;
    return result;
}

} // namespace boundary_dp
