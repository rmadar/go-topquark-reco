package tbuilder

import (
	"math"
	"sort"

	"go-hep.org/x/hep/fmom"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/spatial/r3"
)

// ellipsis returns the t,tbar pair solution according to the Ellipsis method.
func ellipsis(
	lep, lepbar, jet, jetbar fmom.PxPyPzE,
	etx, ety float64,
	debug bool) (r3.Vec, r3.Vec, bool) {

	bldr := newEllipsisBuilder(lep, lepbar, jet, jetbar, etx, ety, debug)
	return bldr.run()
}

type ellipsisBuilder struct {
	lep     fmom.PxPyPzE
	lepbar  fmom.PxPyPzE
	bjet    fmom.PxPyPzE
	bjetbar fmom.PxPyPzE

	etx float64
	ety float64

	debug bool
}

func newEllipsisBuilder(
	lep, lepbar, jet, jetbar fmom.PxPyPzE,
	etx, ety float64,
	//mTop, mTopbar, mw, mwbar, mnu, mnubar float64,
	debug bool) *ellipsisBuilder {
	return &ellipsisBuilder{
		lep:     lep,
		lepbar:  lepbar,
		bjet:    jet,
		bjetbar: jetbar,

		etx: etx,
		ety: ety,

		debug: debug,
	}
}

func (bldr *ellipsisBuilder) run() (r3.Vec, r3.Vec, bool) {

	ps := bldr.neutrinoMomenta()
	if len(ps) == 0 {
		return r3.Vec{}, r3.Vec{}, false
	}

	t, tbar := bldr.selectSolution(ps[0], ps[1])
	return t, tbar, true
}

var (
	zero3x3 = mat.NewDense(3, 3, nil)
)

func (bldr *ellipsisBuilder) neutrinoMomenta() [][]r3.Vec {
	nuEllCalc := newNeutrinoEllipseCalculator(
		bldr.bjet, bldr.lepbar,	mTop, mW, mNu,
	)
	nPerp := nuEllCalc.getNeutrinoEllipse()
	if mat.Equal(nPerp, zero3x3) {
		return nil
	}

	nubarEllCalc := newNeutrinoEllipseCalculator(
		bldr.bjetbar, bldr.lep, mTopbar, mWbar, mNubar,
	)
	nbarPerp := nubarEllCalc.getNeutrinoEllipse()
	if mat.Equal(nbarPerp, zero3x3) {
		return nil
	}

	gamma := mat.NewDense(3, 3, []float64{
		-1, +0, bldr.etx,
		+0, -1, bldr.ety,
		+0, +0, 1,
	})

	nperp := mat.NewDense(3, 3, nil)
	nperp.Mul(nbarPerp, gamma)
	nperp.Mul(gamma.T(), nperp)

	var (
		nuPerps    = intersectEllEll(nPerp, nperp)
		nubarPerps = make([]r3.Vec, 0, len(nuPerps))
	)

	switch len(nuPerps) {
	case 0:
		var (
			hperp = nuEllCalc.hperp
			xp    = mat.NewDense(3, 3, nil)
		)

		xp.Mul(nperp, hperp)
		xp.Mul(hperp.T(), xp)

		d := mat.NewDense(3, 3, []float64{
			0, -1, 0,
			1, +0, 0,
			0, +0, 0,
		})

		xpd := mat.NewDense(3, 3, nil)
		xpd.Mul(xp, d)

		mp := mat.NewDense(3, 3, nil)
		mp.Add(xpd.T(), xpd)

		u := mat.NewDense(3, 3, []float64{
			1, 0, +0,
			0, 1, +0,
			0, 0, -1,
		})

		sols := intersectEllEll(mp, u)
		if len(sols) == 0 {
			return nil
		}

		type ksolT struct {
			k float64
			v r3.Vec
		}
		var (
			ksols = make([]ksolT, 0, len(sols))
			tX    = mat.NewVecDense(3, nil)
			tmp   = make([]float64, 3)
		)
		for i := range sols {
			sol := mat.NewVecDense(3, []float64{
				sols[i].X, sols[i].Y, sols[i].Z,
			})
			for j := 0; j < 3; j++ {
				col := mat.NewVecDense(3, mat.Col(tmp, j, xp))
				tX.SetVec(j, mat.Dot(sol, col))
			}
			ksols = append(ksols, ksolT{
				k: math.Abs(mat.Dot(tX, sol)),
				v: sols[i],
			})
		}
		sort.Slice(ksols, func(i, j int) bool {
			return ksols[i].k < ksols[j].k
		})
		var v mat.VecDense
		v.MulVec(hperp, vecDenseFrom(ksols[0].v))
		nuPerps = append(nuPerps, r3VecFrom(&v))

		var (
			d2      mat.Dense
			npNup   mat.VecDense
			d2NpNup mat.VecDense
		)

		d2.Mul(d, d)
		npNup.MulVec(nPerp, vecDenseFrom(nuPerps[0]))
		d2NpNup.MulVec(&d2, &npNup)

		linePerp := mat.NewVecDense(3, nil)
		{
			nv := nuPerps[0]
			d2 := d2NpNup.RawVector().Data
			linePerp.SetVec(0, nv.Y*d2[2]-nv.Z*d2[1])
			linePerp.SetVec(1, nv.Z*d2[0]-nv.X*d2[2])
			linePerp.SetVec(2, nv.X*d2[1]-nv.Y*d2[0])
		}

		var (
			nuPerpPTmp = intersectEllLine(nil, nperp, r3VecFrom(linePerp))

			dist = 1e9
			idst = -999
		)

		if len(nuPerpPTmp) == 0 {
			return nil
		}

		for i, v := range nuPerpPTmp {
			d := math.Hypot(
				v.X-nuPerps[0].X,
				v.Y-nuPerps[0].Y,
			)
			if d < dist {
				dist = d
				idst = i
			}
		}
		if idst < 0 {
			return nil
		}

		var nubarPerpTmp mat.VecDense
		nubarPerpTmp.MulVec(gamma, vecDenseFrom(nuPerpPTmp[idst]))
		nubarPerps = append(nubarPerps, r3VecFrom(&nubarPerpTmp))

	default:
		vec := mat.NewVecDense(3, nil)
		for _, v := range nuPerps {
			vec.SetVec(0, v.X)
			vec.SetVec(1, v.Y)
			vec.SetVec(2, v.Z)
			vec.MulVec(gamma, vec)
			nubarPerps = append(nubarPerps, r3.Vec{
				X: vec.AtVec(0),
				Y: vec.AtVec(1),
				Z: vec.AtVec(2),
			})
		}
	}

	var (
		hNu           = nuEllCalc.h
		hNubar        = nubarEllCalc.h
		hPerpInvNu    = nuEllCalc.hperpInv
		hPerpInvNubar = nubarEllCalc.hperpInv

		sNu    mat.Dense
		sNubar mat.Dense

		pNus    []r3.Vec
		pNubars []r3.Vec
	)

	sNu.Mul(hNu, hPerpInvNu)
	sNubar.Mul(hNubar, hPerpInvNubar)

	for _, v := range nuPerps {
		var (
			vec = mat.NewVecDense(3, []float64{v.X, v.Y, v.Z})
			nu  mat.VecDense
		)
		nu.MulVec(&sNu, vec)
		pNus = append(pNus, r3.Vec{
			X: nu.AtVec(0),
			Y: nu.AtVec(1),
			Z: nu.AtVec(2),
		})
	}

	for _, v := range nubarPerps {
		var (
			vec = mat.NewVecDense(3, []float64{v.X, v.Y, v.Z})
			nu  mat.VecDense
		)
		nu.MulVec(&sNubar, vec)
		pNubars = append(pNubars, r3.Vec{
			X: nu.AtVec(0),
			Y: nu.AtVec(1),
			Z: nu.AtVec(2),
		})
	}

	return [][]r3.Vec{pNus, pNubars}
}

func (bldr *ellipsisBuilder) selectSolution(nu, nubar []r3.Vec) (top, topbar r3.Vec) {
	mttMin := 100000.0
	if len(nu) != len(nubar) {
		panic("nu/nubar slices len differ")
	}

	for i := range nu {
		v1 := nu[i]
		p1 := fmom.NewPxPyPzE(v1.X, v1.Y, v1.Z, r3.Norm(v1))

		v2 := nubar[i]
		p2 := fmom.NewPxPyPzE(v2.X, v2.Y, v2.Z, r3.Norm(v2))

		var (
			topCalc    fmom.PxPyPzE
			topbarCalc fmom.PxPyPzE
			ttbar      fmom.PxPyPzE
		)

		fmom.IAdd(&topCalc, &p1)
		fmom.IAdd(&topCalc, &bldr.lepbar)
		fmom.IAdd(&topCalc, &bldr.bjet)

		fmom.IAdd(&topbarCalc, &p2)
		fmom.IAdd(&topbarCalc, &bldr.lep)
		fmom.IAdd(&topbarCalc, &bldr.bjetbar)

		fmom.IAdd(&ttbar, &topCalc)
		fmom.IAdd(&ttbar, &topbarCalc)

		if m := ttbar.M(); m < mttMin {
			mttMin = m
			top = fmom.VecOf(&topCalc)
			topbar = fmom.VecOf(&topbarCalc)
		}
	}

	return top, topbar
}

func intersectEllEll(a, b *mat.Dense) []r3.Vec {
	var (
		detA float64
		detB float64
	)

	for j := 0; j < 3; j++ {
		detA += a.At(0, j) * cofactor(a, 0, j)
		detB += b.At(0, j) * cofactor(b, 0, j)
	}

	if math.Abs(detB) > math.Abs(detA) {
		a, b = b, a
	}

	var aInvB mat.Dense
	aInvB.Inverse(a)
	aInvB.Mul(&aInvB, b)

	var eigen mat.Eigen
	if !eigen.Factorize(&aInvB, mat.EigenNone) {
		panic("could not factorize eigen values")
	}
	var (
		iEig  = -1
		eigVs = make([]complex128, 3)
	)
	eigen.Values(eigVs)
	for i, v := range eigVs {
		if imag(v) != 0 {
			continue
		}
		iEig = i
	}

	var m mat.Dense
	m.Scale(-real(eigVs[iEig]), a)
	m.Add(b, &m)

	var (
		lines = factorDegenerate(&m)
		pts   []r3.Vec
	)

	for _, line := range lines {
		pts = intersectEllLine(pts, a, line)
	}

	return pts
}

func intersectEllLine(pts []r3.Vec, a *mat.Dense, line r3.Vec) []r3.Vec {

	var (
		E   = a
		L   = mat.NewVecDense(3, []float64{line.X, line.Y, line.Z})
		LxE = mat.NewDense(3, 3, nil)
		eig mat.Eigen
	)
	// cross product LxE
	for j := 0; j < 3; j++ {
		LxE.Set(0, j, L.AtVec(1)*E.At(2, j)-L.AtVec(2)*E.At(1, j))
		LxE.Set(1, j, L.AtVec(2)*E.At(0, j)-L.AtVec(0)*E.At(2, j))
		LxE.Set(2, j, L.AtVec(0)*E.At(1, j)-L.AtVec(1)*E.At(0, j))
	}

	if !eig.Factorize(LxE, mat.EigenRight) {
		panic("could not factorize eigen vectors")
	}

	type ksolT struct {
		k float64
		v r3.Vec
	}

	var (
		eigs mat.CDense
		sols = make([]ksolT, 0, 3)
	)
	eig.VectorsTo(&eigs)

	for j := 0; j < 3; j++ {
		v := r3.Vec{
			X: real(eigs.At(0, j)),
			Y: real(eigs.At(1, j)),
			Z: real(eigs.At(2, j)),
		}
		v = v.Scale(1 / r3.Norm2(v))
		lv := line.Dot(v)
		var ve r3.Vec
		for jj := 0; jj < 3; jj++ {
			col := r3.Vec{
				X: E.At(0, jj),
				Y: E.At(1, jj),
				Z: E.At(2, jj),
			}
			vv := v.Dot(col)

			switch jj {
			case 0:
				ve.X = vv
			case 1:
				ve.Y = vv
			case 2:
				ve.Z = vv
			}
		}
		vev := ve.Dot(v)
		k := lv*lv + vev*vev
		v = v.Scale(1 / v.Z)
		sols = append(sols, ksolT{
			k: k,
			v: v,
		})
	}
	if len(sols) < 2 {
		return pts
	}

	sort.Slice(sols, func(i, j int) bool {
		return sols[i].k < sols[j].k
	})

	for i := 0; i < 2; i++ {
		sol := sols[i]
		if sol.k > 1e-12 {
			continue
		}
		pts = append(pts, sol.v)
	}

	return pts
}

func cofactor(a *mat.Dense, row, col int) float64 {
	var (
		di, dj int
		m      [2][2]float64
	)
	for i := 0; i < 2; i++ {
		if i == row {
			di = 1
		}
		dj = 0
		for j := 0; j < 2; j++ {
			if j == col {
				dj = 1
			}
			m[i][j] = a.At(i+di, j+dj)
		}
	}
	return ipow(-1, row+col) * (m[0][0]*m[1][1] - m[1][0]*m[0][1])
}

func factorDegenerate(g *mat.Dense) []r3.Vec {
	var (
		lp r3.Vec
		lm r3.Vec
	)

	if g.At(0, 0) == 0 && g.At(1, 1) == 0 {
		lp.X = g.At(0, 1)
		lp.Y = 0
		lp.Z = g.At(1, 2)
		lm.X = 0
		lm.Y = g.At(0, 1)
		lm.Z = g.At(0, 2) - g.At(1, 2)
		return []r3.Vec{lp, lm}
	}

	swapXY := math.Abs(g.At(0, 0)) > math.Abs(g.At(1, 1))
	if swapXY {
		idx := func(i, j int) int {
			return 3*j + i
		}
		sli := g.RawMatrix().Data
		sli[idx(0, 1)], sli[idx(1, 0)] = sli[idx(1, 0)], sli[idx(0, 1)]
		sli[idx(0, 0)], sli[idx(1, 1)] = sli[idx(1, 1)], sli[idx(0, 0)]
		sli[idx(0, 2)], sli[idx(1, 2)] = sli[idx(1, 2)], sli[idx(0, 2)]
		sli[idx(2, 0)], sli[idx(2, 1)] = sli[idx(2, 1)], sli[idx(2, 0)]
	}

	g22 := cofactor(g, 2, 2)

	switch {
	case g22 == 0 && g.At(1, 1) != 0:
		g00 := cofactor(g, 0, 0)
		if g00 > 0 {
			return nil
		}
		s00 := math.Sqrt(-g00)

		lp.X = g.At(0, 1)
		lp.Y = g.At(1, 1)
		lp.Z = g.At(1, 2) + s00
		lm.X = g.At(0, 1)
		lm.Y = g.At(1, 1)
		lm.Z = g.At(1, 2) - s00
	default:
		if g22 > 0 {
			return nil
		}
		var (
			ig22 = 1 / g22
			x0   = cofactor(g, 0, 2) * ig22
			y0   = cofactor(g, 1, 2) * ig22
			s22  = math.Sqrt(-g22)
		)
		lp.X = +g.At(0, 1) + s22
		lp.Y = +g.At(1, 1)
		lp.Z = -g.At(1, 1)*y0 - x0*(g.At(0, 1)+s22)
		lm.X = +g.At(0, 1) - s22
		lm.Y = +g.At(1, 1)
		lm.Z = -g.At(1, 1)*y0 - x0*(g.At(0, 1)-s22)

		if swapXY {
			lp.X, lp.Y = lp.Y, lp.X
			lm.X, lm.Y = lm.Y, lm.X
		}
	}

	return []r3.Vec{lp, lm}
}

func ipow(v float64, i int) float64 {
	if i < 0 {
		panic("ellipsis: negative ipow")
	}
	o := 1.0
	for j := 0; j < i; j++ {
		o *= v
	}
	return o
}

func vecDenseFrom(v r3.Vec) *mat.VecDense {
	return mat.NewVecDense(3, []float64{
		v.X, v.Y, v.Z,
	})
}

func r3VecFrom(v *mat.VecDense) r3.Vec {
	return r3.Vec{
		X: v.AtVec(0),
		Y: v.AtVec(1),
		Z: v.AtVec(2),
	}
}
