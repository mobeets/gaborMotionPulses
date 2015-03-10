function data = loadData(infile)
    load(infile, 'X', 'Y', 'R', 'Xxy');    
    [ny, nt, ns] = size(X);
    Xf = permute(X, [1 3 2]);
    X = reshape(Xf, ny, nt*ns);
    D = asd.sqdist.spaceTime(Xxy, ns, nt);
	Xxy(:,2) = -Xxy(:,2); % flip y-axis
    
    data.Xf = Xf;
    data.X = X;
    data.Y_all = Y;
    data.R = R;
    data.D = D;
    data.ndeltas = size(D, 3);
    data.Xxy = Xxy;
    data.ns = ns;
    data.nt = nt;
end
