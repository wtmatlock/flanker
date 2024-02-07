FROM condaforge/miniforge3:latest
RUN mamba install -c bioconda flanker
WORKDIR /app
COPY . .
RUN pip install pytest
RUN pytest