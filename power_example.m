for i = 1:length(phon_index)-1
   in_F = fft(resample(phon_index(i):phon_index(i+1)));
   out_F = fft(output(phon_index(i):phon_index(i+1)));

   in_pow = in_F.*conj(in_F);
   out_pow = out_F.*conj(out_F);

   total_pow_in = sum(in_pow);
   total_pow_out = sum(out_pow);

   pow_ratio(i) = total_pow_in/total_pow_out;
%     output(phon_index(i):phon_index(i+1)) = output(phon_index(i):phon_index(i+1)) .* pow_ratio(i);

   if isnan(pow_ratio(i)) || isinf(pow_ratio(i)) || pow_ratio(i) > 150 || pow_ratio(i) < 0.3
        output(phon_index(i):phon_index(i+1)) = output(phon_index(i):phon_index(i+1));
   else
       pow_ratio(i);
       output(phon_index(i):phon_index(i+1)) = output(phon_index(i):phon_index(i+1)) .* pow_ratio(i) .* window(@hamming,length(output(phon_index(i):phon_index(i+1))));
%         pow_ratio(i)
%         phon_index(i)
%         phon_index(i+1)
   end
end
