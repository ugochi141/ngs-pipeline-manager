from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from contextlib import asynccontextmanager
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

@asynccontextmanager
async def lifespan(app: FastAPI):
    logger.info("Starting up ngs-pipeline-manager...")
    yield
    logger.info("Shutting down ngs-pipeline-manager...")

app = FastAPI(
    title="ngs-pipeline-manager",
    description="Automated next-generation sequencing data analysis pipeline with bioinformatics workflows",
    version="1.0.0",
    lifespan=lifespan
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.get("/")
async def root():
    return {
        "project": "ngs-pipeline-manager",
        "status": "operational",
        "description": "Automated next-generation sequencing data analysis pipeline with bioinformatics workflows"
    }

@app.get("/health")
async def health_check():
    return {"status": "healthy"}

@app.get("/api/v1/status")
async def api_status():
    return {
        "api_version": "v1",
        "status": "operational",
        "endpoints_available": True
    }
