import React from 'react';

const InputPanel = ({ 
  molA, setMolA, 
  molB, setMolB, 
  config, setConfig, 
  onAnalyze, loading, 
  testCases, onLoadTestCase 
}) => {
  const handleWeightChange = (key, value) => {
    const newAlpha = key === 'alpha' ? parseFloat(value) : config.alpha;
    const newBeta = key === 'beta' ? parseFloat(value) : config.beta;
    
    // Ensure weights sum to 1.0 (simple normalization)
    const total = newAlpha + newBeta;
    setConfig({
      alpha: newAlpha / total,
      beta: newBeta / total,
    });
  };

  return (
    <div className="bg-white rounded-2xl shadow-sm border border-slate-200 p-6 space-y-6">
      <div className="flex items-center gap-2 mb-2">
        <div className="w-1 h-6 bg-accent rounded-full" />
        <h2 className="font-semibold text-slate-700">Molecular Inputs</h2>
      </div>

      <div className="space-y-4">
        <div className="space-y-2">
          <label className="text-xs font-medium text-slate-500 uppercase tracking-wider">Molecule A (SMILES)</label>
          <div className="relative">
            <input 
              type="text"
              value={molA}
              onChange={(e) => setMolA(e.target.value)}
              placeholder="e.g. CC(=O)OC1=CC=CC=C1C(=O)O"
              className="w-full p-3 text-[12px] font-mono bg-slate-50 border border-slate-200 rounded-xl focus:ring-2 focus:ring-accent/20 focus:border-accent outline-none transition-all overflow-x-auto whitespace-nowrap"
            />
            <div className="absolute right-3 bottom-2 text-[10px] text-slate-400 font-mono pointer-events-none">
              {molA.length} chars
            </div>
          </div>
        </div>

        <div className="space-y-2">
          <label className="text-xs font-medium text-slate-500 uppercase tracking-wider">Molecule B / Reference Ligand (SMILES)</label>
          <div className="relative">
            <input 
              type="text"
              value={molB}
              onChange={(e) => setMolB(e.target.value)}
              placeholder="e.g. CC(=O)OC1=CC=CC=C1C(=O)O"
              className="w-full p-3 text-[12px] font-mono bg-slate-50 border border-slate-200 rounded-xl focus:ring-2 focus:ring-accent/20 focus:border-accent outline-none transition-all overflow-x-auto whitespace-nowrap"
            />
            <div className="absolute right-3 bottom-2 text-[10px] text-slate-400 font-mono pointer-events-none">
              {molB.length} chars
            </div>
          </div>
        </div>
      </div>

      <div className="space-y-4 pt-4 border-t border-slate-100">
        <div className="flex items-center gap-2 mb-2">
          <div className="w-1 h-6 bg-accent rounded-full" />
          <h3 className="text-sm font-semibold text-slate-700">Scoring Weights</h3>
        </div>
        
        <div className="space-y-4">
          <div className="space-y-2">
            <div className="flex justify-between text-xs font-medium">
              <span className="text-slate-500">2D Similarity (α)</span>
              <span className="text-accent font-mono">{config.alpha.toFixed(2)}</span>
            </div>
            <input 
              type="range" min="0" max="1" step="0.01"
              value={config.alpha}
              onChange={(e) => handleWeightChange('alpha', e.target.value)}
              className="w-full h-1.5 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-accent"
            />
          </div>

          <div className="space-y-2">
            <div className="flex justify-between text-xs font-medium">
              <span className="text-slate-500">3D Similarity (β)</span>
              <span className="text-accent font-mono">{config.beta.toFixed(2)}</span>
            </div>
            <input 
              type="range" min="0" max="1" step="0.01"
              value={config.beta}
              onChange={(e) => handleWeightChange('beta', e.target.value)}
              className="w-full h-1.5 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-accent"
            />
          </div>
        </div>
      </div>

      <div className="space-y-3 pt-4">
        <button 
          onClick={onAnalyze}
          disabled={loading}
          className={`w-full py-3 px-4 rounded-xl font-semibold text-white transition-all flex items-center justify-center gap-2 ${
            loading ? 'bg-slate-400 cursor-not-allowed' : 'bg-accent hover:bg-blue-700 active:scale-[0.98] shadow-md shadow-accent/20'
          }`}
        >
          {loading ? (
            <>
              <div className="w-4 h-4 border-2 border-white/30 border-t-white rounded-full animate-spin" />
              <span>Analyzing...</span>
            </>
          ) : (
            <>
              <span>Run Analysis</span>
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="M13 10V3L4 14h7v7l9-11h-7z" />
              </svg>
            </>
          )}
        </button>

        <div className="space-y-2">
          <label className="text-xs font-medium text-slate-500 uppercase tracking-wider">Quick Test Cases</label>
          <select 
            onChange={(e) => onLoadTestCase(testCases[e.target.value])}
            defaultValue=""
            className="w-full p-2 text-sm bg-slate-50 border border-slate-200 rounded-lg outline-none focus:ring-2 focus:ring-accent/20"
          >
            <option value="" disabled>Select a case...</option>
            {testCases.map((tc, idx) => (
              <option key={idx} value={idx}>{tc.label}</option>
            ))}
          </select>
        </div>
      </div>
    </div>
  );
};

export default InputPanel;