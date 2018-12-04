FROM rocker/shiny:3.5.0
ENV APP_DIR=/srv/shinyapps/app


# Install required system libraries
RUN sudo apt-get update && apt-get install -y \
    libssl-dev \
    libxml2-dev

# Just copy install.R over to install necessary R packages
COPY app/install.R /tmp/install.R
RUN R -f /tmp/install.R

# Copy Shiny app context to where rocker expects it to reside
COPY app/ $APP_DIR

EXPOSE 3242

WORKDIR $APP_DIR

# Run runShiny.R script
CMD ["sh", "-c", "R -f runShiny.R"]