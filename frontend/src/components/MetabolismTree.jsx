import React, { useEffect, useState, useMemo } from 'react';
import ReactFlow, { Background, Controls, Handle, Position, useNodesState, useEdgesState, MarkerType } from 'reactflow';
import 'reactflow/dist/style.css';
import dagre from 'dagre';
import { Activity, FlaskConical, Info } from 'lucide-react';

// --- Layout Algorithm ---
const getLayoutedElements = (nodes, edges) => {
  const dagreGraph = new dagre.graphlib.Graph();
  dagreGraph.setDefaultEdgeLabel(() => ({}));
  dagreGraph.setGraph({ rankdir: 'TB', ranksep: 280, nodesep: 40 });

  nodes.forEach((node) => dagreGraph.setNode(node.id, { width: 300, height: 260 }));
  edges.forEach((edge) => dagreGraph.setEdge(edge.source, edge.target));

  dagre.layout(dagreGraph);

  return {
    nodes: nodes.map((node) => {
      const pos = dagreGraph.node(node.id);
      return { ...node, position: { x: pos.x - 150, y: pos.y - 130 } };
    }),
    edges
  };
};

// --- Custom Node (Now with Toxicity Colors & Hover Stats) ---
const BioNode = ({ data }) => {
  const renderSVG = (smiles) => {
    if (!window.RDKitModule || !smiles) return null;
    try {
      const mol = window.RDKitModule.get_mol(smiles.trim().replace(/\s/g, ''));
      const svg = mol.get_svg();
      mol.delete();
      return svg;
    } catch (e) { return null; }
  };

  // 1. Determine Shape
  let shapeClass = "rounded-2xl"; 
  if (data.transformation === 'Initial Molecule') shapeClass = "clip-diamond";
  if (data.transformation === 'Phase I' || data.transformation === 'Halogenation') shapeClass = "rounded-full";
  if (data.transformation === 'Phase II' || data.transformation === 'Alkylation') shapeClass = "clip-hexagon";

  // 2. Determine Toxicity Border Color
  const isToxic = data.toxicity >= 0.50;
  const borderColor = isToxic ? "border-red-500" : "border-green-500";
  const shadowColor = isToxic ? "shadow-red-500/30" : "shadow-green-500/30";

  return (
    <div className="flex flex-col items-center group">
      <Handle type="target" position={Position.Top} className="!bg-slate-400 w-4 h-4 border-2 border-white" />
      
      <div 
        onClick={() => data.onClick(data.smiles)}
        className={`relative w-64 h-64 bg-white border-[6px] transition-all duration-300 flex items-center justify-center p-6 cursor-pointer
          ${data.isActive ? `scale-105 z-10 shadow-[0_0_30px] ${shadowColor}` : 'shadow-md hover:scale-105'}
          ${shapeClass} ${borderColor}`}
      >
        <div 
          className="w-full h-full flex items-center justify-center pointer-events-none"
          dangerouslySetInnerHTML={{ __html: renderSVG(data.smiles) || `<span class='text-sm font-bold text-slate-500'>${data.name}</span>` }}
        />
        
        {/* On-Hover Detail Card */}
        <div className="absolute -bottom-20 left-1/2 -translate-x-1/2 opacity-0 group-hover:opacity-100 transition-opacity bg-slate-900 text-white text-xs font-bold py-3 px-5 rounded-2xl whitespace-nowrap z-50 shadow-xl flex flex-col gap-1 items-center">
          <span className={isToxic ? "text-red-400" : "text-green-400"}>
            hERG Risk: {(data.toxicity * 100).toFixed(1)}%
          </span>
          <span className="text-[10px] text-slate-400 font-normal">Click to set as Active Drug</span>
        </div>
      </div>

      <div className="mt-5 text-center bg-white/90 backdrop-blur-sm px-4 py-2 rounded-xl shadow-sm border border-slate-100">
        <p className="text-sm font-black text-slate-900 uppercase tracking-tight">{data.name}</p>
        <p className="text-xs font-bold text-slate-500 mt-1 uppercase tracking-widest">{data.transformation}</p>
      </div>

      <Handle type="source" position={Position.Bottom} className="!bg-slate-400 w-4 h-4 border-2 border-white" />
    </div>
  );
};

// --- Main Component ---
const MetabolismTree = ({ smiles, apiBaseUrl, onNodeClick, activeSmiles }) => {
  const [nodes, setNodes, onNodesChange] = useNodesState([]);
  const [edges, setEdges, onEdgesChange] = useEdgesState([]);
  const [loading, setLoading] = useState(false);
  const [simMode, setSimMode] = useState('liver');
  
  const nodeTypes = useMemo(() => ({ bioNode: BioNode }), []);

  useEffect(() => {
    if (!smiles) return;

    const fetchMetabolism = async () => {
      setLoading(true);
      try {
        const response = await fetch(`${apiBaseUrl}/metabolize`, {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ smiles, max_depth: 2, mode: simMode }),
        });
        const treeData = await response.json();

        const flatNodes = [];
        const flatEdges = [];
        const smilesMap = new Map(); 

        const traverse = (node, parentId = null, edgeLabel = "") => {
          let nodeId = smilesMap.get(node.smiles);
          
          if (!nodeId) {
            nodeId = `node-${smilesMap.size}`;
            smilesMap.set(node.smiles, nodeId);
            flatNodes.push({
              id: nodeId,
              type: 'bioNode',
              data: { 
                smiles: node.smiles, 
                name: node.name, 
                transformation: node.transformation,
                toxicity: node.toxicity, // PASS TOXICITY TO NODE
                onClick: onNodeClick,
                isActive: node.smiles === activeSmiles
              },
            });
          }

          if (parentId) {
            flatEdges.push({
              id: `e-${parentId}-${nodeId}`,
              source: parentId,
              target: nodeId,
              label: edgeLabel,
              labelStyle: { fill: '#1e293b', fontWeight: 800, fontSize: '13px' },
              labelBgPadding: [12, 6],
              labelBgBorderRadius: 6,
              labelBgStyle: { fill: '#ffffff', fillOpacity: 0.9, stroke: '#cbd5e1', strokeWidth: 1 },
              style: { stroke: '#94a3b8', strokeWidth: 2 }, 
              markerEnd: { type: MarkerType.ArrowClosed, color: '#94a3b8', width: 20, height: 20 },
            });
          }

          if (node.children) {
            node.children.forEach(child => traverse(child, nodeId, child.name));
          }
        };

        traverse(treeData);
        const { nodes: lNodes, edges: lEdges } = getLayoutedElements(flatNodes, flatEdges);
        setNodes(lNodes);
        setEdges(lEdges);
      } catch (err) {
        console.error(err);
      } finally {
        setLoading(false);
      }
    };

    fetchMetabolism();
  }, [smiles, apiBaseUrl, onNodeClick, activeSmiles, simMode]);

  return (
    <div className="mt-4 space-y-4">
      <style>{`
        .clip-diamond { clip-path: polygon(50% 0%, 100% 50%, 50% 100%, 0% 50%); }
        .clip-hexagon { clip-path: polygon(25% 0%, 75% 0%, 100% 50%, 75% 100%, 25% 100%, 0% 50%); }
        .react-flow__attribution { display: none; }
      `}</style>
      
      <div className="flex flex-col xl:flex-row items-start xl:items-center justify-between gap-4 bg-white p-4 rounded-2xl border border-slate-200 shadow-sm">
        
        {/* Dynamic Header */}
        <div className="flex items-center gap-3">
          <div className={`p-2 rounded-lg shadow-lg ${simMode === 'liver' ? 'bg-blue-600 shadow-blue-200' : 'bg-purple-600 shadow-purple-200'}`}>
            {simMode === 'liver' ? <Activity size={20} className="text-white" /> : <FlaskConical size={20} className="text-white" />}
          </div>
          <div>
            <h3 className="text-base font-black text-slate-900 uppercase tracking-widest">
              {simMode === 'liver' ? 'Metabolic Trajectory' : 'MedChem Generator'}
            </h3>
            <p className="text-xs text-slate-500 font-bold uppercase">
              {simMode === 'liver' ? 'Hepatic Clearance Simulation' : 'Generative Multi-Parameter Optimization'}
            </p>
          </div>
        </div>

        {/* --- UI LEGEND (Shapes & Colors) --- */}
        <div className="hidden md:flex items-center gap-4 bg-slate-50 px-4 py-2 rounded-xl border border-slate-200 text-[10px] font-bold text-slate-600 uppercase tracking-wider">
           <div className="flex items-center gap-2"><div className="w-3 h-3 bg-slate-300 clip-diamond" /> Parent</div>
           <div className="flex items-center gap-2"><div className="w-3 h-3 bg-slate-300 rounded-full" /> Phase I / Minor</div>
           <div className="flex items-center gap-2"><div className="w-3 h-3 bg-slate-300 clip-hexagon" /> Phase II / Major</div>
           <div className="w-[1px] h-4 bg-slate-300 mx-1" />
           <div className="flex items-center gap-2"><div className="w-3 h-3 border-2 border-green-500 rounded-full" /> Safe (&lt;50%)</div>
           <div className="flex items-center gap-2"><div className="w-3 h-3 border-2 border-red-500 rounded-full" /> Toxic (&ge;50%)</div>
        </div>

        {/* Toggle Switch */}
        <div className="flex bg-slate-100 p-1 rounded-xl border border-slate-200 shrink-0">
          <button 
            onClick={() => setSimMode('liver')}
            className={`px-4 py-2 text-xs font-bold uppercase tracking-widest rounded-lg transition-all ${
              simMode === 'liver' ? 'bg-white text-blue-600 shadow-sm border border-slate-200' : 'text-slate-400 hover:text-slate-600'
            }`}
          >
            Liver Simulator
          </button>
          <button 
            onClick={() => setSimMode('medchem')}
            className={`px-4 py-2 text-xs font-bold uppercase tracking-widest rounded-lg transition-all ${
              simMode === 'medchem' ? 'bg-white text-purple-600 shadow-sm border border-slate-200' : 'text-slate-400 hover:text-slate-600'
            }`}
          >
            MedChem Opt
          </button>
        </div>

      </div>

      <div className="h-[900px] w-full bg-slate-50 rounded-3xl border-2 border-slate-200 shadow-inner relative overflow-hidden">
        {loading && (
          <div className="absolute inset-0 z-50 bg-white/80 backdrop-blur-sm flex flex-col items-center justify-center">
            <div className={`w-16 h-16 border-4 border-t-transparent rounded-full animate-spin mb-6 ${simMode === 'liver' ? 'border-blue-600' : 'border-purple-600'}`} />
            <p className="text-sm font-black text-slate-800 uppercase tracking-widest">Generating Network...</p>
          </div>
        )}
        
        <ReactFlow
          nodes={nodes}
          edges={edges}
          onNodesChange={onNodesChange}
          onEdgesChange={onEdgesChange}
          nodeTypes={nodeTypes}
          fitView
          fitViewOptions={{ padding: 0.1, minZoom: 0.2 }} 
          minZoom={0.1}
          maxZoom={1.5}
          panOnScroll={true}
          zoomOnScroll={false}
          panOnDrag={true}
        >
          <Background variant="dots" gap={30} size={2} color="#cbd5e1" />
          <Controls className="!bg-white !border-slate-200 !shadow-xl !fill-slate-700" showInteractive={false} />
        </ReactFlow>
      </div>
    </div>
  );
};

export default MetabolismTree;