import React from 'react';
import { Radar, RadarChart, PolarGrid, PolarAngleAxis, PolarRadiusAxis, ResponsiveContainer, Legend, Tooltip } from 'recharts';

const DescriptorRadar = ({ results }) => {
  if (!results) return null;

  const { descriptors_a = {}, descriptors_b = {} } = results;

  const normalize = (val, min, max) => {
    if (val === undefined || val === null || isNaN(parseFloat(val))) return 0
    const n = parseFloat(val)
    return Math.max(0, Math.min(100, ((n - min) / (max - min)) * 100))
  }

  const radarData = [
    {
      subject: "LogP",
      A: normalize(descriptors_a?.logP, -2, 8),
      B: normalize(descriptors_b?.logP, -2, 8),
      fullMark: 100
    },
    {
      subject: "MW",
      A: normalize(descriptors_a?.MW, 0, 800),
      B: normalize(descriptors_b?.MW, 0, 800),
      fullMark: 100
    },
    {
      subject: "TPSA",
      A: normalize(descriptors_a?.TPSA, 0, 150),
      B: normalize(descriptors_b?.TPSA, 0, 150),
      fullMark: 100
    },
    {
      subject: "HBD",
      A: normalize(descriptors_a?.HBD, 0, 8),
      B: normalize(descriptors_b?.HBD, 0, 8),
      fullMark: 100
    },
    {
      subject: "HBA",
      A: normalize(descriptors_a?.HBA, 0, 12),
      B: normalize(descriptors_b?.HBA, 0, 12),
      fullMark: 100
    }
  ];

  return (
    <div className="bg-white rounded-xl border border-slate-200 overflow-hidden shadow-sm p-4">
      <div className="h-[300px] w-full">
        <ResponsiveContainer width="100%" height="100%">
          <RadarChart cx="50%" cy="50%" outerRadius="75%" data={radarData}>
            <PolarGrid stroke="#e2e8f0" />
            <PolarAngleAxis dataKey="subject" tick={{ fill: '#64748b', fontSize: 11 }} />
            <PolarRadiusAxis angle={30} domain={[0, 100]} tick={false} axisLine={false} />
            <Radar
              name="Molecule A"
              dataKey="A"
              stroke="#2563eb"
              fill="#2563eb"
              fillOpacity={0.3}
              strokeWidth={2}
            />
            <Radar
              name="Molecule B"
              dataKey="B"
              stroke="#16a34a"
              fill="#16a34a"
              fillOpacity={0.3}
              strokeWidth={2}
            />
            <Tooltip formatter={(val, name) => [val.toFixed(1), name]} />
            <Legend verticalAlign="bottom" height={36}/>
          </RadarChart>
        </ResponsiveContainer>
      </div>
    </div>
  );
};

export default DescriptorRadar;