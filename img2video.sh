#!/bin/bash

# 定义图片和视频的根目录
IMG_ROOT="Video&Img/Image"
VIDEO_ROOT="Video&Img/Video"

# 确保视频根目录存在
mkdir -p "$VIDEO_ROOT"

# 遍历图片根目录下的所有子文件夹
for IMG_FOLDER in "$IMG_ROOT"/*/; do
    # 获取子文件夹的名称
    FOLDER_NAME=$(basename "$IMG_FOLDER")
    
    # 定义输出视频文件的路径
    VIDEO_PATH="$VIDEO_ROOT/$FOLDER_NAME.mp4"
    
    # 使用 ffmpeg 将图片合成视频
    ffmpeg -framerate 24 -pattern_type glob -i "$IMG_FOLDER/image-*.png" -c:v libx264 -pix_fmt yuv420p -b:v 24M "$VIDEO_PATH"
    
    # 输出处理信息
    echo "Processed $IMG_FOLDER and saved video to $VIDEO_PATH"
done

echo "All videos have been processed and saved to $VIDEO_ROOT"