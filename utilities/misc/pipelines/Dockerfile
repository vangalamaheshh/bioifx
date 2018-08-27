FROM mvangala/base-nginx:python3.6

MAINTAINER Mahesh Vangala <vangalamaheshh@gmail.com>

RUN set -ex \
  && apt-get update -y \
  && apt-get install -y sudo \
  && echo "deb http://ppa.launchpad.net/webupd8team/java/ubuntu xenial main" | \
  sudo tee /etc/apt/sources.list.d/webupd8team-java.list \
  && echo "deb-src http://ppa.launchpad.net/webupd8team/java/ubuntu xenial main" | \
  sudo tee -a /etc/apt/sources.list.d/webupd8team-java.list \
  && apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys EEA14886 \
  && apt-get update -y \
  && echo oracle-java8-installer shared/accepted-oracle-license-v1-1 select true | debconf-set-selections \
  && apt-get install -y oracle-java8-installer \
  && apt-get install -y maven

ENV LC_ALL "C.UTF-8"
ENV LANG "C.UTF-8"

# install dockerflow
WORKDIR /usr/local/bin
RUN set -ex \
  && git clone https://github.com/vangalamaheshh/my-beam.git \
  && cd my-beam/variant-caller \
  && mvn package -DskipTests

# install flask and bigquery client api
RUN set -ex \
  && pip install flask 

RUN set -ex \
  && mkdir -p /usr/local/bin/pipelines-api 

WORKDIR /usr/local/bin/pipelines-api

COPY uwsgi.ini ./
COPY entrypoint.sh /
COPY server.py ./

# set up flask with nginx
# By default, allow unlimited file sizes, modify it to limit the file sizes
# To have a maximum of 1 MB (Nginx's default) change the line to:
# ENV NGINX_MAX_UPLOAD 1m
ENV NGINX_MAX_UPLOAD 0

# Which uWSGI .ini file should be used, to make it customizable
ENV UWSGI_INI /usr/local/bin/pipelines-api/uwsgi.ini
# URL under which static (not modified by Python) files will be requested
# They will be served by Nginx directly, without being handled by uWSGI
ENV STATIC_URL /static
# Absolute path in where the static files wil be
ENV STATIC_PATH /usr/local/bin/pipelines-api/static
# If STATIC_INDEX is 1, serve / with /static/index.html directly (or the static URL configured)
# ENV STATIC_INDEX 1
ENV STATIC_INDEX 0
# Copy the entrypoint that will generate Nginx additional configs
COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]

CMD ["/usr/bin/supervisord"]

