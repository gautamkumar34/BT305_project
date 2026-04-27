import React, { useEffect, useRef, useState } from 'react';

const StructureViewer = ({ molA, molB, alignmentData, loading }) => {
  const viewerRef = useRef(null);
  const viewerInstance = useRef(null);
  const [rdkit, setRdkit] = useState(null);
  const [svgA, setSvgA] = useState(null);
  const [svgB, setSvgB] = useState(null);

  // Initialize RDKit.js
  useEffect(() => {
    const loadRDKit = async () => {
      if (window.RDKitModule) {
        setRdkit(window.RDKitModule);
        return;
      }
      if (window.initRDKitPromise) {
        await window.initRDKitPromise;
        setRdkit(window.RDKitModule);
      }
    };
    loadRDKit();
  }, []);

  // Update 2D SVGs when SMILES change
  useEffect(() => {
    if (!rdkit) return;

    const renderSVG = (smiles) => {
      if (!smiles) return null;
      const cleanSmiles = smiles.trim().replace(/\s/g, '');
      const mol = rdkit.get_mol(cleanSmiles);
      if (!mol || !mol.is_valid()) return null;
      const svg = mol.get_svg();
      mol.delete();
      return svg;
    };

    setSvgA(renderSVG(molA));
    setSvgB(renderSVG(molB));
  }, [rdkit, molA, molB]);

  // Update 3D Viewer when alignment data changes
  useEffect(() => {
    if (!alignmentData || !viewerRef.current) return;

    try {
      if (!viewerInstance.current) {
        viewerInstance.current = window.$3Dmol.createViewer(viewerRef.current, {
          backgroundColor: 'white',
        });
      }

      const viewer = viewerInstance.current;
      viewer.clear();

      // Molecule A - Blue
      viewer.addModel(alignmentData.mol_a_sdf, 'sdf');
      viewer.setStyle({ model: 0 }, { stick: { radius: 0.15, color: '#2563eb' } });

      // Molecule B - Green
      viewer.addModel(alignmentData.mol_b_sdf, 'sdf');
      viewer.setStyle({ model: 1 }, { stick: { radius: 0.15, color: '#16a34a' } });

      viewer.zoomTo();
      viewer.render();
    } catch (err) {
      console.error('3Dmol rendering error:', err);
    }
  }, [alignmentData]);

  return (
    <div className="space-y-6">
      {/* 2D View Section - Side by Side Boxes */}
      <div className="grid grid-cols-2 gap-4">
        <div className="flex flex-col items-center p-4 bg-white rounded-xl border border-slate-100 shadow-sm relative overflow-hidden">
          <div className="absolute top-0 left-0 w-1 h-full bg-blue-500" />
          <span className="text-[10px] font-bold text-slate-400 uppercase tracking-widest mb-4">DRUG A (REFERENCE)</span>
          <div 
            className="w-full h-40 flex items-center justify-center p-2"
            dangerouslySetInnerHTML={{ __html: svgA || '<div class="text-slate-300 text-[10px]">NO STRUCTURE</div>' }}
          />
        </div>
        <div className="flex flex-col items-center p-4 bg-white rounded-xl border border-slate-100 shadow-sm relative overflow-hidden">
          <div className="absolute top-0 left-0 w-1 h-full bg-green-500" />
          <span className="text-[10px] font-bold text-slate-400 uppercase tracking-widest mb-4">DRUG B (ANALOG)</span>
          <div 
            className="w-full h-40 flex items-center justify-center p-2"
            dangerouslySetInnerHTML={{ __html: svgB || '<div class="text-slate-300 text-[10px]">NO STRUCTURE</div>' }}
          />
        </div>
      </div>

      {/* 3D View Section */}
      <div className="bg-white p-4 rounded-xl border border-slate-100 shadow-sm">
        <div className="flex justify-between items-center mb-3">
          <h3 className="text-xs font-bold text-slate-700 uppercase tracking-wider">3D Flexible Alignment</h3>
          <div className="flex gap-4">
            <div className="flex items-center gap-1.5 text-[9px] font-bold text-slate-500 uppercase tracking-tighter">
              <div className="w-2 h-2 rounded-full bg-blue-600" /> Drug A
            </div>
            <div className="flex items-center gap-1.5 text-[9px] font-bold text-slate-500 uppercase tracking-tighter">
              <div className="w-2 h-2 rounded-full bg-green-600" /> Drug B
            </div>
          </div>
        </div>
        
        <div className="relative group">
          <div 
            ref={viewerRef} 
            className="w-full h-[400px] rounded-lg bg-slate-50 border border-slate-100 overflow-hidden"
            style={{ position: 'relative' }}
          />
          
          {!alignmentData && !loading && (
            <div className="absolute inset-0 flex items-center justify-center bg-slate-50/80 backdrop-blur-[1px]">
              <p className="text-slate-400 text-xs font-semibold uppercase tracking-widest">Awaiting Analysis...</p>
            </div>
          )}
          
          {loading && (
            <div className="absolute inset-0 flex items-center justify-center bg-white/60 backdrop-blur-sm z-20">
              <div className="flex flex-col items-center gap-2">
                <div className="w-8 h-8 border-4 border-blue-600 border-t-transparent rounded-full animate-spin" />
                <span className="text-[10px] font-bold text-slate-600 uppercase tracking-wider">Performing Alignment...</span>
              </div>
            </div>
          )}
        </div>
      </div>
    </div>
  );
};

export default StructureViewer;