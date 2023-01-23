To run this file:

docker build -t my-jupyter-notebook .
docker run -p 8888:8888 my-jupyter-notebook

---- R studio
docker run -it -d -p 7777:8787 --name tviblindi_container -v ~/M2_CLIP_visit:/data -e PASSWORD=pass1234 stuchly/tviblindi
