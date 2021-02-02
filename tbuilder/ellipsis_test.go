package tbuilder

import (
	"testing"

	"go-hep.org/x/hep/fmom"
	"gonum.org/v1/gonum/spatial/r3"
)

func TestEllipsis(t *testing.T) {
	for _, tc := range []struct {
		lep1, lep2 fmom.PxPyPzE
		jet1, jet2 fmom.PxPyPzE
		etx, ety   float64

		top1, top2 r3.Vec
		want       bool
	}{
		{
			lep1: fmom.NewPxPyPzE(-12.608633, -31.456136, -155.817925, 159.460663),
			lep2: fmom.NewPxPyPzE(-60.243233, -43.445790, -78.368434, 107.974101),
			jet1: fmom.NewPxPyPzE(43.576636, 4.794894, -101.073176, 110.467331),
			jet2: fmom.NewPxPyPzE(50.678302, 66.619242, -125.304938, 151.095413),
			etx:  17.232637,
			ety:  1.596588,

			want: true,
			top1: r3.Vec{-47.7501392810703, -8.384996581583039, -194.27801797016332},
			top2: r3.Vec{86.38584828107031, 6.493794581583032, -347.5347243998728},
		},
		{
			lep1: fmom.NewPxPyPzE(-56.843705, -49.294411, 168.804950, 184.814114),
			lep2: fmom.NewPxPyPzE(-12.150408, 51.595663, 49.461753, 72.499723),
			jet1: fmom.NewPxPyPzE(-44.565618, -50.361997, 283.036653, 290.955627),
			jet2: fmom.NewPxPyPzE(24.923737, -8.670943, 27.339602, 38.411682),
			etx:  7.855825738005946e+01,
			ety:  7.536658086618009e+01,

			want: false,
		},
		{
			lep1: fmom.NewPxPyPzE(-56.843705, -49.294411, 168.804950, 184.814114),
			lep2: fmom.NewPxPyPzE(-12.150408, 51.595663, 49.461753, 72.499723),
			jet1: fmom.NewPxPyPzE(24.923737, -8.670943, 27.339602, 38.411682),
			jet2: fmom.NewPxPyPzE(-44.565618, -50.361997, 283.036653, 290.955627),
			etx:  7.855825738005946e+01,
			ety:  7.536658086618009e+01,

			want: true,
			top1: r3.Vec{139.07785507872507, 378.0718815484892, 402.2881339969903},
			top2: r3.Vec{-220.45703800579722, -325.98621951617133, 795.6366258567978},
		},
	} {
		t.Run("", func(t *testing.T) {
			t1, t2, ok := ellipsis(
				tc.lep1, tc.lep2, tc.jet1, tc.jet2,
				tc.etx, tc.ety,
				false,
			)

			if got, want := ok, tc.want; got != want {
				t.Fatalf("invalid reco decision got=%v, want=%v", got, want)
			}

			if got, want := t1, tc.top1; got != want {
				t.Fatalf("invalid top1 3-vector: got=%v, want=%v", got, want)
			}

			if got, want := t2, tc.top2; got != want {
				t.Fatalf("invalid top2 3-vector: got=%v, want=%v", got, want)
			}
		})
	}
}
