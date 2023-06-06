# Builds celligner docker image
# Note that this docker image does not have taigapy or mnnpy installed

docker buildx build --platform linux/amd64 --push -t us.gcr.io/bdekegel/celligner:latest .