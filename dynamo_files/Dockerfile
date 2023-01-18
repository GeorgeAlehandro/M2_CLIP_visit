FROM jupyter/minimal-notebook:latest
# Set the app directory as the working directory
WORKDIR /app

# Add the current directory to the app directory
ADD . /app

RUN pip install -r requirements.txt
