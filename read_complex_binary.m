function v = read_complex_binary (filename, count)

  % usage: read_complex_binary (filename, [count])
  %
  %  open filename and return the contents as a column vector,
  %  treating them as 32 bit complex numbers
  %

  m = nargchk (1,2,nargin);
  if (m)
    usage (m);
  end

  if (nargin < 2)
    count = Inf;
  end

  f = fopen (filename, 'rb');
  if (f < 0)
    v = 0;
  else
    t = fread (f, [2, count], 'float');
    fclose (f);
    v = t(1,:) + t(2,:)*i;
    [r, c] = size (v);
    v = reshape (v, c, r);
  end
end