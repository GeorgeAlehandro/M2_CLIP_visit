# VeloSIM data generation
## Use on R>4.1

docker build -t r-version-4-2-2 .

docker run --rm -p 8888:8787 -v /home/georgealehandro/M2_CLIP_visit:/data -e PASSWORD=password r-version-4-2-2
