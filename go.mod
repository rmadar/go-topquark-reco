module github.com/rmadar/go-topquark-reco

go 1.14

require (
	github.com/rmadar/go-lorentz-vector v0.0.0-20200408140142-c838293f9a81
	go-hep.org/x/exp v0.7.1
	go-hep.org/x/hep v0.28.4-0.20201127142150-e07f188986ac
	golang.org/x/exp v0.0.0-20200513190911-00229845015e
	gonum.org/v1/gonum v0.8.2
)

replace gonum.org/v1/gonum => github.com/gonum/gonum v0.8.1-0.20201202015544-18580834651d
