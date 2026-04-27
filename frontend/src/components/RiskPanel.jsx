import React from 'react';
import { AlertTriangle, CheckCircle, ShieldAlert } from 'lucide-react';
import { RISK_LEVELS } from '../constants';

const RiskPanel = ({ riskA, rulesA, riskB, rulesB, mlToxicity }) => {
  const renderDrugCard = (label, risk, rules, prob) => {
    let badgeColor = "bg-green-500";
    if (risk === "HIGH_RISK") badgeColor = "bg-red-500";
    if (risk === "MODERATE_RISK") badgeColor = "bg-yellow-500";
    
    return (
      <div className="bg-white p-4 rounded-xl border border-slate-200 flex flex-col h-full shadow-sm">
        <div className="flex justify-between items-center mb-2">
          <span className="text-[10px] font-bold text-slate-400 uppercase tracking-widest">{label}</span>
          <div className={`px-2 py-0.5 rounded-full text-[9px] font-bold uppercase text-white ${badgeColor}`}>
            {risk || 'LOW_RISK'}
          </div>
        </div>

        {prob !== undefined && !isNaN(parseFloat(prob)) && (
          <div className="mb-4 text-[11px] font-bold text-slate-700">
            {(parseFloat(prob) * 100).toFixed(1)}% RISK PROB
          </div>
        )}
        
        <div className="flex-grow">
          <span className="text-[10px] font-bold text-slate-500 uppercase block mb-2 tracking-tight">Risk Factors</span>
          <div className="space-y-2">
            {rules && rules.length > 0 ? (
              rules.map((rule, idx) => (
                <div key={idx} className="flex items-start gap-2 text-xs text-slate-600 leading-relaxed font-medium">
                  <div className={`w-1 h-1 rounded-full mt-1.5 shrink-0 ${badgeColor}`} />
                  <span>{rule}</span>
                </div>
              ))
            ) : (
              <span className="text-xs text-slate-400 italic">No specific risk flags detected</span>
            )}
          </div>
        </div>
      </div>
    );
  };

  return (
    <div className="space-y-4">
      <div className="px-4 py-3 border-b border-slate-100 bg-slate-50/50">
        <h3 className="text-xs font-bold text-slate-500 uppercase tracking-wider">Toxicity Risk Assessment</h3>
      </div>
      
      <div className="grid grid-cols-2 gap-4">
        {renderDrugCard('Drug A', riskA, rulesA, mlToxicity?.prob_a)}
        {renderDrugCard('Drug B', riskB, rulesB, mlToxicity?.prob_b)}
      </div>
      
      <div className="text-center">
        <p className="text-[10px] text-slate-400 italic">
          Based on dual-similarity scoring and known cardiotoxicity patterns.
        </p>
      </div>
    </div>
  );
};

export default RiskPanel;