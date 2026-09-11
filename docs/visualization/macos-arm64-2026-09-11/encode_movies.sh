#!/bin/sh
# Encode the 4D voxel frames to MP4. Recorded so the movies are reproducible.
# usage: encode_movies.sh <anim_dir> <out_dir>
set -e
A="$1"; O="$2"
for seq in parcels quorum; do
  ffmpeg -y -framerate 12 -i "$A/anim_${seq}_%03d.png" \
    -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" \
    -c:v libx264 -pix_fmt yuv420p -crf 20 "$O/voxels_4d_${seq}.mp4"
done
