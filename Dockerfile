# Build stage: compile C++ extension
FROM python:3.11-slim AS builder

# Install build dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy source code
COPY . .

# Install Python package (non-editable mode)
RUN pip install pybind11 numpy && \
    pip install ./python/

# Runtime stage: smaller image
FROM python:3.11-slim

WORKDIR /app

# Copy built extension from builder
COPY --from=builder /usr/local/lib/python3.11/site-packages /usr/local/lib/python3.11/site-packages

# Copy app files only (not the whole build)
COPY app/ /app/app/

# Install runtime dependencies
RUN pip install streamlit matplotlib numpy

# Expose Streamlit port
EXPOSE 8501

# Run Streamlit
ENTRYPOINT ["streamlit", "run", "app/streamlit_app.py", "--server.port=8501", "--server.address=0.0.0.0"]