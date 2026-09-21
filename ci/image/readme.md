Build:
```shell
docker build -t geosets-image .
```
or
```shell
docker build -t geosets-image-python -f ci/image/Dockerfile.python .
```

Run:
```shell
docker run -it geosets-image /bin/bash
```
