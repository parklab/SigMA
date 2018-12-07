FROM rocker/shiny:3.5.0
ENV APP_DIR=/srv/shinyapps/app


# Install required system libraries
RUN sudo apt-get update && apt-get install -y \
    libssl-dev \
    libxml2-dev

# Just copy install.R over to install necessary R packages
COPY ./install.R /tmp/install.R
RUN R -f /tmp/install.R

WORKDIR $APP_DIR
EXPOSE 3242

# Copy SigMa app context to where rocker expects it to reside
COPY . $APP_DIR

# Run shiny app
CMD ["R", "-e", "shiny::runApp('/srv/shinyapps/app/shiny/', 3242, host='0.0.0.0')"]
