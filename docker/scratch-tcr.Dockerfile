# Use a specific version of Ubuntu as the base image
FROM --platform=linux/x86_64 rocker/verse:latest

# Set the working directory inside the container
WORKDIR /opt

# Timezone settings
ENV TZ=US/Central
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
    echo $TZ > /etc/timezone

# Install system dependencies
RUN apt-get update && apt-get install -y \
    software-properties-common \
    dirmngr \
    gnupg \
    apt-transport-https \
    ca-certificates \
    wget \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev

# Updating quarto to Quarto v1.4.553
RUN wget https://github.com/quarto-dev/quarto-cli/releases/download/v1.4.553/quarto-1.4.553-linux-amd64.deb -O quarto-1.4.553-linux-amd64.deb
RUN dpkg -i quarto-1.4.553-linux-amd64.deb

# Install Python3
RUN apt-get install -y \
    python3 \
    python3-pip

# Install Python packages for data science
RUN python3 -m pip install --no-cache-dir numpy pandas scikit-learn matplotlib seaborn jupyter
RUN python3 -m pip install --no-cache-dir jupyter-cache
RUN python3 -m pip install --no-cache-dir papermill

RUN python3 -m pip install --no-cache-dir scanpy scirpy session_info

# Cleaning apt-get cache
RUN apt-get clean
RUN rm -rf /var/lib/apt/lists/*

# Command to run on container start
CMD ["cellranger"]
