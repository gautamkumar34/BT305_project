FROM python:3.11-slim

# Install system dependencies for RDKit and Math libraries
RUN apt-get update && apt-get install -y \
    libxrender1 libxext6 libfontconfig1 \
    libopenblas-dev \
    && rm -rf /var/lib/apt/lists/*

# Create a non-root user (Required by Hugging Face)
RUN useradd -m -u 1000 user
USER user
ENV PATH="/home/user/.local/bin:$PATH"
WORKDIR /app

# Set ML Environment Variables
ENV OPENBLAS_NUM_THREADS=1 \
    OMP_NUM_THREADS=1 \
    TF_USE_LEGACY_KERAS=1

# Copy and install dependencies
COPY --chown=user requirements.txt .
RUN pip install --no-cache-dir --upgrade pip
RUN pip install --no-cache-dir -r requirements.txt

# Copy the rest of your clean code
COPY --chown=user . .

# Hugging Face Spaces must listen on port 7860
EXPOSE 7860

# Start the app (Point to your main.py and the 'app' object)
CMD ["uvicorn", "src.api:app", "--host", "0.0.0.0", "--port", "7860"]