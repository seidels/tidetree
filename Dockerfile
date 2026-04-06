# Dockerfile to build container for unit testing
FROM debian:stable
RUN apt-get update
RUN apt-get install -y openjdk-21-jdk openjfx ant git

WORKDIR /root
RUN git clone https://github.com/CompEvol/beast2.git ../beast2
RUN git clone https://github.com/CompEvol/BeastFX.git ../BeastFX
RUN git clone https://github.com/tgvaughan/feast.git ../feast
ADD . ./
ENTRYPOINT JAVA_FX_HOME=/usr/share/java/ ant test	