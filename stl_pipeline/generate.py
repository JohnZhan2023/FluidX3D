import trimesh
import numpy as np
import os

# 输出文件夹
output_path = "generated_stl/"
os.makedirs(output_path, exist_ok=True)

# 生成形状的函数
def create_shapes():
    shapes = {}

    # **1. 立方体变体**
    for i, scale in enumerate([(1,1,1), (1,2,1), (2,1,1), (1,1,2)]):
        mesh = trimesh.creation.box(extents=scale)
        mesh.apply_transform(trimesh.transformations.rotation_matrix(np.pi/4, [1,0,0]))  # 旋转
        shapes[f"cube_{i}"] = mesh

    # **2. 球体变体**
    for i, sub in enumerate([2, 3, 4, 5]):
        shapes[f"sphere_{i}"] = trimesh.creation.icosphere(subdivisions=sub, radius=1)

    # **3. 圆柱体变体**
    for i, (r, h) in enumerate([(0.5, 2), (1, 1.5), (0.7, 2.5), (1, 3)]):
        mesh = trimesh.creation.cylinder(radius=r, height=h)
        mesh.apply_transform(trimesh.transformations.rotation_matrix(np.pi/6 * i, [0,1,0]))
        shapes[f"cylinder_{i}"] = mesh

    # **4. 长方体**
    for i, dims in enumerate([(1, 2, 3), (2, 3, 1), (3, 1, 2)]):
        shapes[f"rectangular_{i}"] = trimesh.creation.box(extents=dims)

    # **5. 椭球体**
    for i, scales in enumerate([(1, 1.5, 1), (1.2, 1, 0.8), (1, 2, 1)]):
        mesh = trimesh.creation.icosphere(subdivisions=3, radius=1)
        mesh.apply_scale(scales)
        shapes[f"ellipsoid_{i}"] = mesh

    # **6. 锥体**
    for i, h in enumerate([2, 3, 4, 5]):
        shapes[f"cone_{i}"] = trimesh.creation.cone(radius=1, height=h)

    # # **7. 截锥体**
    # for i, (r1, r2, h) in enumerate([(1, 0.5, 2), (1, 0.3, 3), (0.8, 0.4, 2.5)]):
    #     shapes[f"frustum_{i}"] = trimesh.creation.conical_frustum(radius_top=r2, radius_base=r1, height=h)


    # **8. 环面（甜甜圈）**
    for i, r in enumerate([(1, 0.3), (1.2, 0.4), (1.5, 0.5)]):
        mesh = trimesh.creation.torus(r[0], r[1])
        mesh.apply_transform(trimesh.transformations.rotation_matrix(np.pi/3 * i, [1,0,0]))
        shapes[f"torus_{i}"] = mesh

    return shapes

# 生成 40 个 STL 文件
shapes = create_shapes()
for name, shape in shapes.items():
    shape.export(f"{output_path}{name}.stl")
    print(f"已生成 {output_path}{name}.stl")
