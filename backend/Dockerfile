FROM python:3.11-slim

# System libs required by scientific wheels (scanpy/leidenalg/Pillow/hdf5).
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential libhdf5-dev libgl1 && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY app ./app

EXPOSE 8000
# Visium uploads can be large; uvicorn has no body cap of its own, but front
# a reverse proxy with a matching client_max_body_size in production.
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000"]
