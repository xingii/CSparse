# CloudCompare 运行版（直接对应 Pan GitHub 的 cloudFilter 主流程）

你要求的是“只把 GitHub 代码改成 CloudCompare 可运行”。

本插件现在只保留并迁移了 `RoadMarkingExtraction` 中最核心、最独立的 `cloudFilter` 思路：

- OTSU 强度阈值分割
- SOR 统计离群点剔除

即：**OTSU + SOR**（与 Pan 仓库 `pointcloudprocess::cloudFilter` 一致的主干逻辑）。

---

## 1) 编译

```bash
cmake -S CloudCompareRoadMarkingPlugin \
  -B build-roadmarking \
  -DCloudCompare_DIR=/your/cloudcompare/install/lib/cmake/CloudCompare

cmake --build build-roadmarking -j
```

---

## 2) 运行

1. 把动态库放到 CloudCompare 的 `plugins/` 目录。
2. 打开 CloudCompare，加载点云。
3. 选中点云并把 intensity 设为当前显示标量场。
4. 菜单执行：`Plugins > Extract Road Markings (Pan 2019)`。
5. 输出 `*_road_markings`。

---

## 3) 代码映射

- `PanRoadMarkingAdapter::extractMask`：对应 OTSU 阈值分割。
- `PanRoadMarkingAdapter::sorFilter`：对应 SOR 去噪。
- `qRoadMarking::doAction`：CloudCompare 插件入口（读取当前选中点云 + 结果回写 DB Tree）。

