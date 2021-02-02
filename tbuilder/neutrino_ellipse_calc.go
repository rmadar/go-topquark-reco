package tbuilder

import (
	"fmt"
	"math"

	"go-hep.org/x/hep/fmom"
	"gonum.org/v1/gonum/mat"
)

type neutrinoEllipseCalculator struct {
	bjet fmom.PxPyPzE
	lep  fmom.PxPyPzE

	bjetBeta, bjetBeta2, bjetGamma, bjetGamma2 float64
	lepBeta, lepBeta2, lepGamma, lepGamma2     float64

	// particle masses
	mw, mt, mnu    float64
	mw2, mt2, mnu2 float64

	// numbers
	x0, x0p  float64
	sx, sy   float64
	epsilon2 float64

	cos, sin float64 //cosine and sine of theta_{b,mu}

	omega  float64
	Omega  float64
	x1, y1 float64
	z2     float64

	ab *mat.Dense
	al *mat.Dense

	ht       *mat.Dense
	h        *mat.Dense
	hperp    *mat.Dense
	hperpInv *mat.Dense

	nperp *mat.Dense
}

func newNeutrinoEllipseCalculator(
	bjet, lepbar fmom.PxPyPzE,
	mt, mw, mnu float64,
) *neutrinoEllipseCalculator {
	nec := &neutrinoEllipseCalculator{
		bjet:      bjet,
		lep:       lepbar,
		bjetBeta:  betaOf(&bjet),
		bjetGamma: gammaOf(&bjet),
		lepBeta:   betaOf(&lepbar),
		lepGamma:  gammaOf(&lepbar),
		mw:        mw,
		mt:        mt,
		mnu:       mnu,
		mw2:       mw * mw,
		mt2:       mt * mt,
		mnu2:      mnu * mnu,
		cos:       fmom.CosTheta(&lepbar, &bjet),
		ab:        mat.NewDense(4, 4, nil),
		al:        mat.NewDense(4, 4, nil),
		ht:        mat.NewDense(3, 3, nil),
		h:         mat.NewDense(3, 3, nil),
		hperp:     mat.NewDense(3, 3, nil),
		hperpInv:  mat.NewDense(3, 3, nil),
		nperp:     mat.NewDense(3, 3, nil),
	}
	nec.bjetBeta2 = nec.bjetBeta * nec.bjetBeta
	nec.bjetGamma2 = nec.bjetGamma * nec.bjetGamma
	nec.lepBeta2 = nec.lepBeta * nec.lepBeta
	nec.lepGamma2 = nec.lepGamma * nec.lepGamma

	nec.sin = math.Sqrt(1 - nec.cos*nec.cos)
	return nec
}

func (nec *neutrinoEllipseCalculator) getNeutrinoEllipse() *mat.Dense {
	nec.Wsurface()
	nec.leptonEllipsoid()
	nec.bJetEllipsoid()
	nec.neutrinoSolution()
	nec.labSystemTransform()
	return nec.nperp // FIXME(sbinet): clone?
}

func (nec *neutrinoEllipseCalculator) Wsurface() {
	nec.x0p = -(0.5 / nec.bjet.E()) * (nec.mt2 - nec.mw2 - nec.bjet.M2())
	nec.x0 = -(0.5 / nec.lep.E()) * (nec.mw2 - nec.lep.M2() - nec.mnu2)
	nec.sx = (1. / nec.lepBeta2) * (nec.x0*nec.lepBeta - nec.lep.P()*(1.-nec.lepBeta2))
	nec.epsilon2 = (1. - nec.lepBeta2) * (nec.mw2 - nec.mnu2)
}

func (nec *neutrinoEllipseCalculator) leptonEllipsoid() {
	nec.al.Set(0, 0, 1.-nec.lepBeta2)
	nec.al.Set(1, 0, 0)
	nec.al.Set(2, 0, 0)
	nec.al.Set(3, 0, nec.sx*nec.lepBeta2)

	nec.al.Set(0, 1, 0)
	nec.al.Set(1, 1, 1)
	nec.al.Set(2, 1, 0)
	nec.al.Set(3, 1, 0)

	nec.al.Set(0, 2, 0)
	nec.al.Set(1, 2, 0)
	nec.al.Set(2, 2, 1)
	nec.al.Set(3, 2, 0)

	nec.al.Set(0, 3, nec.sx*nec.lepBeta2)
	nec.al.Set(1, 3, 0)
	nec.al.Set(2, 3, 0)
	nec.al.Set(3, 3, nec.mw2-nec.x0*nec.x0-nec.epsilon2)
}

func (nec *neutrinoEllipseCalculator) bJetEllipsoid() {
	nec.ab.Set(0, 0, 1-nec.cos*nec.cos*nec.bjetBeta2)
	nec.ab.Set(1, 0, -nec.cos*nec.sin*nec.bjetBeta2)
	nec.ab.Set(2, 0, 0)
	nec.ab.Set(3, 0, nec.cos*nec.x0p*nec.bjetBeta)

	nec.ab.Set(0, 1, -nec.cos*nec.sin*nec.bjetBeta2)
	nec.ab.Set(1, 1, 1-nec.sin*nec.sin*nec.bjetBeta2)
	nec.ab.Set(2, 1, 0)
	nec.ab.Set(3, 1, nec.sin*nec.x0p*nec.bjetBeta)

	nec.ab.Set(0, 2, 0)
	nec.ab.Set(1, 2, 0)
	nec.ab.Set(2, 2, 1)
	nec.ab.Set(3, 2, 0)

	nec.ab.Set(0, 3, nec.cos*nec.x0p*nec.bjetBeta)
	nec.ab.Set(1, 3, nec.sin*nec.x0p*nec.bjetBeta)
	nec.ab.Set(2, 3, 0)
	nec.ab.Set(3, 3, nec.mw2-nec.x0p*nec.x0p)
}

func (nec *neutrinoEllipseCalculator) neutrinoSolution() {
	nec.sy = (1. / nec.sin) * (nec.x0p/nec.bjetBeta - nec.cos*nec.sx)
	nec.omega = (1. / nec.sin) * (nec.lepBeta/nec.bjetBeta - nec.cos) //only the positive slope
	nec.Omega = math.Sqrt(math.Max(0, nec.omega*nec.omega+1.-nec.lepBeta2))
	Omega2 := nec.Omega * nec.Omega
	nec.x1 = nec.sx - (1./Omega2)*(nec.sx+nec.omega*nec.sy)
	nec.y1 = nec.sy - (1./Omega2)*nec.omega*(nec.sx+nec.omega*nec.sy)
	nec.z2 = nec.x1*nec.x1*Omega2 - (nec.sy-nec.omega*nec.sx)*(nec.sy-nec.omega*nec.sx) - (nec.mw2 - nec.x0*nec.x0 - nec.epsilon2)
	Z := math.Sqrt(math.Max(0, nec.z2))

	nec.ht.Set(0, 0, Z/nec.Omega)
	nec.ht.Set(0, 1, 0)
	nec.ht.Set(0, 2, nec.x1-nec.lep.P())

	nec.ht.Set(1, 0, nec.omega*Z/nec.Omega)
	nec.ht.Set(1, 1, 0)
	nec.ht.Set(1, 2, nec.y1)

	nec.ht.Set(2, 0, 0)
	nec.ht.Set(2, 1, Z)
	nec.ht.Set(2, 2, 0)
}

func (nec *neutrinoEllipseCalculator) labSystemTransform() {
	Rz := rotationMatrix(2, -nec.lep.Phi())
	Ry := rotationMatrix(1, 0.5*math.Pi-thetaOf(&nec.lep))
	bJetP := fmom.VecOf(&nec.bjet)
	bJet_xyz := mat.NewVecDense(3, []float64{bJetP.X, bJetP.Y, bJetP.Z})

	var rM mat.Dense
	rM.Mul(Rz, bJet_xyz)
	rM.Mul(Ry, &rM)

	rA := rM.RawMatrix().Data
	phi := -math.Atan2(rA[2], rA[1])

	Rx := rotationMatrix(0, phi)

	var R = mat.NewDense(3, 3, nil)
	R.Mul(Ry.T(), Rx.T())
	R.Mul(Rz.T(), R)

	nec.h.Reset()
	nec.h.Mul(R, nec.ht)

	h := nec.h.RawMatrix().Data
	hp := nec.hperp.RawMatrix().Data

	copy(hp[:6], h[:6])
	hp[6] = 0
	hp[7] = 0
	hp[8] = 1

	det := mat.Det(nec.hperp) // FIXME(sbinet): work in log-space?
	if det == 0 {
		return
	}

	copy(nec.hperpInv.RawMatrix().Data, hp)
	nec.hperpInv.Inverse(nec.hperpInv)

	U := mat.NewDense(3, 3, []float64{
		1, 0, 0,
		0, 1, 0,
		0, 0, -1,
	})

	nec.nperp.Reset()
	nec.nperp.Mul(U, nec.hperpInv)
	nec.nperp.Mul(nec.hperpInv.T(), nec.nperp)

	//	Rz := r3.NewRotation(-nec.lepton_.Phi(), r3.Vec{Z: 1})
	//	Ry := r3.NewRotation(0.5*math.Pi-thetaOf(&nec.lepton_), r3.Vec{Y: 1})
	//	bJetP := fmom.VecOf(&nec.bJet_)
	//	//bJet_xyz := mat.NewVecDense(3, []float64{bJetP.X, bJetP.Y, bJetP.Z})
	//	rA := Ry.Rotate(Rz.Rotate(bJetP))
	//	phi := -math.Atan2(rA.Z, rA.Y)
	//
	//	Rz := r3.NewRotation(phi, r3.Vec{X:1})
	//	R := Rz.Rotate(Ry.Rotate(
	//	var rM mat.Dense
	//	rM.Mul(Rz, bJet_xyz)
}

func betaOf(p4 fmom.P4) float64 {
	return p4.P() / p4.E()
}

func gammaOf(p4 fmom.P4) float64 {
	b := betaOf(p4)
	return 1.0 / math.Sqrt(1-b*b)
}

func thetaOf(p4 fmom.P4) float64 {
	var (
		px = p4.Px()
		py = p4.Py()
		pz = p4.Pz()
	)
	if px == 0 && py == 0 && pz == 0 {
		return 0
	}
	return math.Atan2(p4.Pt(), pz)
}

func rotationMatrix(axis int, angle float64) *mat.Dense {
	sin, cos := math.Sincos(angle)
	switch axis {
	case 0:
		return mat.NewDense(3, 3, []float64{
			1, 0, 0,
			0, +cos, -sin,
			0, +sin, +cos,
		})
	case 1:
		return mat.NewDense(3, 3, []float64{
			+cos, 0, +sin,
			0, 1, 0,
			-sin, 0, +cos,
		})
	case 2:
		return mat.NewDense(3, 3, []float64{
			+cos, -sin, 0,
			+sin, +cos, 0,
			0, 0, 1,
		})
	default:
		panic(fmt.Errorf("invalid axis=%d", axis))
	}
}
