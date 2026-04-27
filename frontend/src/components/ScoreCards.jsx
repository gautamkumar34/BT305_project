import React from 'react';

const ScoreCard = ({ label, value, isPrimary = false }) => {
  const percentage = Math.min(100, Math.max(0, value * 100));
  
  let textColor = "text-blue-600";
  let barColor = "bg-blue-600";
  
  if (value < 0.4) {
    textColor = "text-red-500";
    barColor = "bg-red-500";
  } else if (value <= 0.7) {
    textColor = "text-yellow-600";
    barColor = "bg-yellow-500";
  }

  if (isPrimary) {
    return (
      <div className="bg-[#0f172a] text-white p-6 rounded-2xl shadow-xl border border-slate-700 relative overflow-hidden">
        <div className="relative z-10">
          <div className="text-[10px] font-bold text-slate-400 uppercase tracking-widest mb-1">{label}</div>
          <div className="text-4xl font-bold font-mono tracking-tight">{value.toFixed(3)}</div>
        </div>
        <div className="absolute top-0 right-0 w-24 h-24 bg-blue-500/10 rounded-full -mr-12 -mt-12" />
      </div>
    );
  }
  
  return (
    <div className="bg-white p-4 rounded-xl border border-slate-100 shadow-sm">
      <div className="flex justify-between items-end mb-2">
        <span className="text-[10px] font-bold text-slate-400 uppercase tracking-widest leading-none">
          {label}
        </span>
        <span className={`text-2xl font-bold font-mono leading-none ${textColor}`}>
          {value.toFixed(3)}
        </span>
      </div>
      <div className="w-full h-1.5 bg-slate-100 rounded-full overflow-hidden">
        <div 
          className={`h-full transition-all duration-1000 ease-out ${barColor}`}
          style={{ width: `${percentage}%` }}
        />
      </div>
    </div>
  );
};

const ScoreCards = ({ tanimoto2D, shape3D, finalScore, o3a_score }) => {
  return (
    <div className="space-y-6">
      <div className="flex items-center gap-2 px-1">
        <div className="w-1 h-5 bg-blue-600 rounded-full" />
        <h2 className="text-xs font-bold text-slate-700 uppercase tracking-wider">Similarity Metrics</h2>
      </div>
      
      <div className="space-y-4">
        <ScoreCard label="TANIMOTO 2D" value={tanimoto2D} />
        <ScoreCard label="SHAPE 3D" value={shape3D} />
        <ScoreCard label="COMBINED SCORE" value={finalScore} isPrimary={true} />
      </div>

      <div className="px-1">
        <div className="text-slate-400 text-xs font-medium">
          O3A Score: {o3a_score?.toFixed(4)}
        </div>
        <div className="text-slate-400 text-[10px] mt-0.5 italic">
          Higher = better alignment
        </div>
      </div>
    </div>
  );
};

export default ScoreCards;