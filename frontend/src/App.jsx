import React, { useState, useEffect } from 'react';
import { analyzeSimilarity, checkHealth } from './api';
import { DEFAULT_CONFIG, RISK_LEVELS, TEST_CASES } from './constants';

// Sub-components
import InputPanel from './components/InputPanel';
import StructureViewer from './components/StructureViewer';
import ScoreCards from './components/ScoreCards';
import RiskPanel from './components/RiskPanel';
import ExplanationBox from './components/ExplanationBox';
import DescriptorTable from './components/DescriptorTable';
import DescriptorRadar from './components/DescriptorRadar';
import MetabolismTree from './components/MetabolismTree'; 

const App = () => {
  const [molA, setMolA] = useState('');
  const [molB, setMolB] = useState('');
  const [config, setConfig] = useState(DEFAULT_CONFIG);
  const [results, setResults] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [backendStatus, setBackendStatus] = useState('checking');

  // Safe fallback for API URL
  const apiBaseUrl = import.meta.env.VITE_API_BASE_URL || 'https://gautam34-biotransformer-api.hf.space';

  useEffect(() => {
    const verifyBackend = async () => {
      const status = await checkHealth();
      setBackendStatus(status.status === 'online' ? 'online' : 'offline');
    };
    verifyBackend();
  }, []);

  const handleAnalyze = async () => {
    if (!molA || !molB) {
      setError('Please provide both SMILES strings');
      return;
    }

    setLoading(true);
    setError(null);
    try {
      const data = await analyzeSimilarity({
        smiles_a: molA,
        smiles_b: molB,
        alpha: config.alpha,
        beta: config.beta,
      });
      setResults(data);
    } catch (err) {
      setError(err instanceof Error ? err.message : String(err));
    } finally {
      setLoading(false);
    }
  };

  const loadTestCase = (testCase) => {
    setMolA(testCase.smiles_a);
    setMolB(testCase.smiles_b);
    setError(null);
  };

  return (
    <div className="min-h-screen bg-slate-50 font-sans text-slate-900">
      {/* Header */}
      <header className="bg-navy text-white py-4 px-6 shadow-lg flex justify-between items-center z-50 relative">
        <div className="flex items-center gap-3">
          <div className="w-10 h-10 bg-accent rounded-lg flex items-center justify-center font-bold text-xl">BT</div>
          <h1 className="text-xl font-bold tracking-tight">
            BioTransformer <span className="font-normal text-slate-300 text-sm ml-2">| Molecular Similarity & Cardiotoxicity</span>
          </h1>
        </div>
        <div className="flex items-center gap-2">
          <div className={`w-3 h-3 rounded-full ${backendStatus === 'online' ? 'bg-green-400' : 'bg-red-400'} animate-pulse`} />
          <span className="text-xs font-mono uppercase tracking-widest text-slate-300">
            Backend: {backendStatus}
          </span>
        </div>
      </header>

      {/* Main Layout: 75% Left/Center (Scrolling) | 25% Right (Sticky) */}
      <main className="p-6 max-w-[1800px] mx-auto flex flex-col xl:flex-row gap-6 items-start">
        
        {/* ========================================== */}
        {/* LEFT & CENTER COLUMN (75% Width) - Scrolls */}
        {/* ========================================== */}
        <div className="w-full xl:w-3/4 flex flex-col gap-6">
          
          {/* 1. Inputs */}
          <InputPanel 
            molA={molA} setMolA={setMolA} 
            molB={molB} setMolB={setMolB} 
            config={config} setConfig={setConfig}
            onAnalyze={handleAnalyze}
            loading={loading}
            testCases={TEST_CASES}
            onLoadTestCase={loadTestCase}
          />
          {error && (
            <div className="p-4 bg-red-50 border border-red-200 text-red-700 rounded-xl text-sm animate-in fade-in slide-in-from-top-2">
              <strong>Error:</strong> {error}
            </div>
          )}

          {/* 2. Molecular Alignment Viewer */}
          <div className="bg-white rounded-2xl shadow-sm border border-slate-200 overflow-hidden flex flex-col">
            <div className="px-6 py-4 border-b border-slate-100 bg-slate-50/50 flex justify-between items-center">
              <h2 className="text-xs font-bold text-slate-700 uppercase tracking-wider">Molecular Alignment</h2>
              <span className="text-[10px] text-slate-400 font-mono">RDKit.js + 3Dmol.js</span>
            </div>
            <div className="p-6">
              <StructureViewer 
                molA={molA} 
                molB={molB} 
                alignmentData={results?.alignment} 
                loading={loading} 
              />
            </div>
          </div>

          {/* 3. Similarity Metrics (Placed under alignment per your sketch) */}
          {results && (
            <ScoreCards 
              tanimoto2D={results.tanimoto_2d} 
              shape3D={results.shape_3d} 
              finalScore={results.final_score} 
              o3a_score={results.o3a_score}
            />
          )}

          {/* 4. Metabolic Trajectory (Full 75% Width) */}
          <div className="w-full">
            {molA ? (
              <MetabolismTree 
                smiles={molA} 
                activeSmiles={molA}
                apiBaseUrl={apiBaseUrl} 
                onNodeClick={(metaboliteSmiles) => {
                  setMolA(metaboliteSmiles);
                  handleAnalyze();
                  window.scrollTo({ top: 0, behavior: 'smooth' });
                }}
              />
            ) : (
              <div className="mt-2 p-12 bg-white rounded-2xl shadow-sm border border-dashed border-slate-200 text-center">
                <p className="text-sm font-bold text-slate-400 uppercase tracking-widest">
                  Enter Drug A to simulate Metabolic Trajectory
                </p>
              </div>
            )}
          </div>
        </div>

        {/* ========================================== */}
        {/* RIGHT SIDEBAR (25% Width) - Sticky */}
        {/* ========================================== */}
        <div className="w-full xl:w-1/4 flex flex-col gap-6 sticky top-6">
          {results ? (
            <>
              {/* Toxicity Assessment */}
              <RiskPanel 
                riskA={results.risk_flag_a} 
                rulesA={results.risk_rules_a} 
                riskB={results.risk_flag_b} 
                rulesB={results.risk_rules_b} 
                mlToxicity={results.ml_toxicity}
              />
              
              {/* AI Interpretation */}
              <ExplanationBox explanation={results.explanation} />
              
              {/* Physicochemical Profile */}
              <div className="bg-white p-6 rounded-2xl shadow-sm border border-slate-200">
                <div className="flex items-center gap-2 px-1 mb-4">
                  <div className="w-1 h-5 bg-blue-600 rounded-full" />
                  <h3 className="text-sm font-semibold text-slate-700">Physicochemical Profile</h3>
                </div>
                <DescriptorTable results={results} />
              </div>

              {/* Property Fingerprint Radar */}
              <div className="bg-white p-6 rounded-2xl shadow-sm border border-slate-200">
                <div className="flex items-center gap-2 px-1 mb-4">
                  <div className="w-1 h-5 bg-blue-600 rounded-full" />
                  <h3 className="text-sm font-semibold text-slate-700">Property Fingerprint</h3>
                </div>
                <DescriptorRadar results={results} />
              </div>
            </>
          ) : (
            <div className="h-[600px] flex flex-col items-center justify-center p-8 text-center bg-white rounded-2xl border-2 border-dashed border-slate-200 text-slate-400">
              <div className="w-16 h-16 mb-4 rounded-full bg-slate-100 flex items-center justify-center text-2xl">📊</div>
              <p className="text-sm font-bold">Analysis Sidebar</p>
              <p className="text-xs mt-2">Run an analysis to view Toxicity, AI Interpretation, and Physicochemical Properties.</p>
            </div>
          )}
        </div>

      </main>
    </div>
  );
};

export default App;