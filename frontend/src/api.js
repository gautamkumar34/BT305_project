import axios from 'axios';
import { API_BASE_URL, ENDPOINTS } from './constants';

const api = axios.create({
  baseURL: API_BASE_URL,
  headers: {
    'Content-Type': 'application/json',
  },
});

/**
 * Perform ligand-based similarity analysis
 * @param {Object} params - Analysis parameters
 * @param {string} params.smiles_a - SMILES of molecule A
 * @param {string} params.smiles_b - SMILES of molecule B
 * @param {number} params.alpha - Weight for 2D similarity
 * @param {number} params.beta - Weight for 3D similarity
 */
export const analyzeSimilarity = async ({ smiles_a, smiles_b, alpha, beta }) => {
  try {
    const response = await api.post(ENDPOINTS.SIMILARITY, {
      smiles_a,
      smiles_b,
      alpha,
      beta,
    });
    return response.data;
  } catch (error) {
    console.error('Error analyzing similarity:', error);
    throw error.response?.data?.detail || new Error('Failed to analyze molecular similarity');
  }
};

/**
 * Validate against known toxicity cases
 * @param {string} case_id - Identifier for the validation case
 */
export const validateToxicity = async (case_id) => {
  try {
    const response = await api.get(`${ENDPOINTS.VALIDATE}/${case_id}`);
    return response.data;
  } catch (error) {
    console.error('Error validating toxicity:', error);
    throw error.response?.data?.detail || new Error('Failed to validate toxicity case');
  }
};

/**
 * Check backend health status
 */
export const checkHealth = async () => {
  try {
    const response = await api.get(ENDPOINTS.HEALTH);
    return response.data;
  } catch (error) {
    console.error('Backend health check failed:', error);
    return { status: 'offline' };
  }
};

export const predictMetabolism = async ({ smiles, max_depth, mode }) => {
  try {
    const response = await api.post(ENDPOINTS.METABOLIZE, {
      smiles,
      max_depth,
      mode,
    });
    return response.data;
  } catch (error) {
    console.error('Error predicting metabolism:', error);
    throw error.response?.data?.detail || new Error('Failed to generate metabolic network');
  }
};

export default api;