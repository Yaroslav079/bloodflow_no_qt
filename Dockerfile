FROM ubuntu:22.04 AS build

RUN apt-get update && \
    apt-get install -y cmake \
                libeigen3-dev \
                libomp-dev \
                build-essential

WORKDIR /bloodflow

COPY ext/ ./ext/
COPY lib/ ./lib/
COPY src/ ./src/
COPY CMakeLists.txt .

RUN cmake -DCMAKE_BUILD_TYPE=Release CMakeLists.txt && \
    cmake --build . --parallel 8

FROM ubuntu:22.04

WORKDIR /app

COPY data/ ./data/

RUN apt-get update && \
    apt-get install -y libomp-dev \
                libgomp1

COPY --from=build \
    ./bloodflow/bloodflow \
    .

# ENTRYPOINT [ "./bloodflow", "-manual", "25", "120", "80", "80", "60" ]

# these defaults should be overwritten by passing to run command list of parameters
# TODO : will it work with sif container?
CMD [ "./bloodflow", "-manual", "25", "120", "80", "80", "60" ]