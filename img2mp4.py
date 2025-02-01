import os
import ffmpeg

# 设置输入文件夹和输出文件夹
input_dir = "C:\\Users\\86177\\Desktop\\workspace\\FluidX3D\\SynData\\test\\images\\0, 0, 0; Zorro character"
output_dir = ""

# 确保输出目录存在
os.makedirs(output_dir, exist_ok=True)

# 遍历输入目录中的每个子文件夹
for subfolder in os.listdir(input_dir):
    subfolder_path = os.path.join(input_dir, subfolder)
    
    # 只处理子文件夹
    if os.path.isdir(subfolder_path):
        # 定义输出视频文件路径
        output_video = os.path.join(output_dir, f"{subfolder}.mp4")
        
        # 使用 ffmpeg-python 创建视频
        input_pattern = os.path.join(subfolder_path, "image-*.png")  # 根据文件名模式选择所有图片
        ffmpeg.input(input_pattern, pattern_type='glob', framerate=24).output(output_video, vcodec='libx264', pix_fmt='yuv420p').run()
        
        print(f"视频已保存为 {output_video}")
