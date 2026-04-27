import React from 'react';

const DescriptorTable = ({ results }) => {
  if (!results) return null;

  const { descriptors_a = {}, descriptors_b = {}, descriptor_delta = {} } = results;

  const safeNum = (val, decimals=2) => {
    const n = parseFloat(val)
    return isNaN(n) ? null : n.toFixed(decimals)
  }

  const formatDelta = (val, type) => {
    if (val === undefined || val === null) return <span className="text-slate-300">—</span>;
    
    const num = parseFloat(val);
    if (isNaN(num)) return <span className="text-slate-300">—</span>;

    const formatted = num > 0 ? `+${num.toFixed(type === 'fraction_csp3' ? 3 : 2)}` : num.toFixed(type === 'fraction_csp3' ? 3 : 2);
    
    let colorClass = 'text-slate-400'; // neutral
    if (type === 'tpsa') {
      colorClass = num > 0 ? 'text-green-600' : 'text-red-500';
    } else if (type === 'logp') {
      colorClass = num < 0 ? 'text-green-600' : 'text-red-500';
    } else if (type === 'hbd') {
      colorClass = num > 0 ? 'text-orange-500' : 'text-slate-400';
    }

    return <span className={colorClass}>{formatted}</span>;
  };

  const rows = [
    {
      label: "Mol Weight (Da)",
      a: safeNum(descriptors_a?.MW),
      b: safeNum(descriptors_b?.MW),
      delta: descriptor_delta?.mw,
      deltaType: "neutral"
    },
    {
      label: "LogP",
      a: safeNum(descriptors_a?.logP),
      b: safeNum(descriptors_b?.logP),
      delta: descriptor_delta?.logp,
      deltaType: "logp"
    },
    {
      label: "H-Bond Donors",
      a: safeNum(descriptors_a?.HBD),
      b: safeNum(descriptors_b?.HBD),
      delta: descriptor_delta?.hbd,
      deltaType: "hbd"
    },
    {
      label: "H-Bond Acceptors",
      a: safeNum(descriptors_a?.HBA),
      b: safeNum(descriptors_b?.HBA),
      delta: descriptor_delta?.hba,
      deltaType: "neutral"
    },
    {
      label: "TPSA (Å²)",
      a: safeNum(descriptors_a?.TPSA, 1),
      b: safeNum(descriptors_b?.TPSA, 1),
      delta: descriptor_delta?.tpsa,
      deltaType: "tpsa"
    },
    {
      label: "Rotatable Bonds",
      a: null,
      b: null,
      delta: descriptor_delta?.rotatable_bonds,
      deltaType: "neutral"
    },
    {
      label: "Fraction CSP3",
      a: null,
      b: null,
      delta: safeNum(descriptor_delta?.fraction_csp3, 3),
      deltaType: "neutral"
    }
  ];

  return (
    <div className="bg-white rounded-xl border border-slate-200 overflow-hidden shadow-sm">
      <div className="overflow-x-auto">
        <table className="w-full text-left text-sm">
          <thead>
            <tr className="bg-slate-50 text-slate-400 font-medium border-b border-slate-100">
              <th className="px-4 py-3 font-semibold">Property</th>
              <th className="px-4 py-3 font-semibold text-right">Mol A</th>
              <th className="px-4 py-3 font-semibold text-right">Mol B</th>
              <th className="px-4 py-3 font-semibold text-right">Δ</th>
            </tr>
          </thead>
          <tbody className="divide-y divide-slate-100">
            {rows.map((row, idx) => (
              <tr key={idx} className="hover:bg-slate-50 transition-colors">
                <td className="px-4 py-3 font-medium text-slate-600">{row.label}</td>
                <td className="px-4 py-3 text-right font-mono text-slate-500">
                  {row.a ?? '—'}
                </td>
                <td className="px-4 py-3 text-right font-mono text-slate-500">
                  {row.b ?? '—'}
                </td>
                <td className="px-4 py-3 text-right font-mono">
                  {formatDelta(row.delta, row.deltaType)}
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
};

export default DescriptorTable;