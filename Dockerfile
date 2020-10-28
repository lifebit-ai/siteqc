FROM nfcore/base@sha256:2043dea2e3215a32576e2e9fa957d8d41f439d209abbf1f858fd02829d2b7d64
#FROM nfcore/base:1.10.2

LABEL authors="Christina Chatzipantsiou, Vlad Dembrovskyi" \
      description="Docker image containing all software requirements for the siteqc pipeline"

# Install GAWK and other needed tools
RUN apt-get update && \
    apt-get install -y \
                   gawk \
                   wget \
                   zip \
                   procps && \
    rm -rf /var/lib/apt/lists/*


# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/siteqc-1.0dev/bin:$PATH

# Install stringi R package and the ones that depend on it.
# (Issue with stringi package from conda that it depends on libicu64 that
# is not available for Debian 10)

RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/stringi/stringi_1.4.5.tar.gz', repos=NULL, type='source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/stringr/stringr_1.3.1.tar.gz', repos=NULL, type='source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/plyr/plyr_1.8.5.tar.gz', repos=NULL, type='source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/reshape2/reshape2_1.4.3.tar.gz', repos=NULL, type='source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/tidyr/tidyr_1.0.2.tar.gz', repos=NULL, type='source')"

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-siteqc-1.0dev > nf-core-siteqc-1.0dev.yml

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron

# Install plink2 Alpha 2.3 Developmental (20 Oct 2020)
# Correct version of Plink2 is not available on Conda
RUN wget http://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20201020.zip && \
    unzip plink2_linux_x86_64_20201020.zip -d plink2 && \
    rm plink2_linux_x86_64_20201020.zip

ENV PATH /plink2:$PATH
