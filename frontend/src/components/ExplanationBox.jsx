import React from 'react';

const ExplanationBox = ({ explanation }) => {
  return (
    <div className="bg-white rounded-2xl shadow-sm border border-slate-200 p-6 space-y-3">
      <div className="flex items-center gap-2 mb-2">
        <div className="w-1 h-6 bg-accent rounded-full" />
        <h2 className="font-semibold text-slate-700">AI Interpretation</h2>
      </div>
      <div className="p-4 bg-slate-50 rounded-xl border border-slate-100">
        <p className="text-sm text-slate-600 leading-relaxed italic">
          "{explanation}"
        </p>
      </div>
      <div className="flex items-center gap-2 text-[10px] text-slate-400 font-medium uppercase tracking-wider">
        <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="M13 16h-1v-4h-1m1-4h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z" />
        </svg>
        Analysis based on structural overlap and pharmacophore mapping
      </div>
    </div>
  );
};

export default ExplanationBox;