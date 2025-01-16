#!/bin/bash

# 设置目标目录和输出文件
VIDEO_DIR="Val_Video&Img/videos"
OUTPUT_FILE="Val_Video&Img/videos.txt"
PROMPT_FILE="Val_Video&Img/prompts.txt"

# 检查目标目录是否存在
if [ ! -d "$VIDEO_DIR" ]; then
    echo "目录 $VIDEO_DIR 不存在，请检查路径。"
    exit 1
fi

# 清空或创建 video.txt 文件
> "$OUTPUT_FILE"

# 遍历目录下的视频文件并写入相对路径到 video.txt
find "$VIDEO_DIR" -type f \( -name "*.mp4" -o -name "*.avi" -o -name "*.mkv" -o -name "*.mov" \) | while read -r file; do
    # 转换为相对路径
    relative_path="${file#$VIDEO_DIR/}"
    echo "videos/$relative_path" >> "$OUTPUT_FILE"
    echo "Visualize the Q-criterion response of the object in the given image to represent its turbulent flow field in fluid dynamics." >> "$PROMPT_FILE"  # 写入空行

done

echo "视频文件路径已写入 $OUTPUT_FILE"
