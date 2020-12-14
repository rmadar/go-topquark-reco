package sonn

// Options encodes the various settings to pass to
// the reconstruction.
type Option func(cfg *config)

// config contains all the possible options and their values.
type config struct {

	// Seed for random numbers
	rndseed uint64  

	// Number of smearing iteration
	smearN int
	
	// Select smeared quantities
	smearAll      bool
	smearLepPt    bool
	smearLepTheta bool
	smearLepAzimu bool
	smearJetPt    bool
	smearJetTheta bool
	smearJetAzimu bool

	// Debug printing
	debug bool
}

// newConfig returns a config type with a set of default options.
func newConfig() *config {
	cfg := &config{
		rndseed:       1234,
		smearAll:      true,
		smearLepPt:    true,
		smearLepTheta: true,
		smearLepAzimu: true,
		smearJetPt:    true,
		smearJetTheta: true,
		smearJetAzimu: true,
		smearN:        10,
		debug:         false,
	}
	return cfg
}

// WithRndSeed sets the seed for random number
// generations used in the smearing. Default is 1234
func WithRndSeed(s uint64) Option {
	return func(cfg *config) {
		cfg.rndseed = s
	}
}

// WithSmearAll enables smearing of all kinematics quantities
// This include lepton pT, polar, and azimuth angles as well
// as jet energy, polar and azimuth angles. Enabled by default.
// Default number of iterations is 10, but can be changed
// using WithSmearN(n) function.
func WithSmearAll(d bool) Option {
	return func(cfg *config) {
		cfg.smearAll = d
	}
}

// WithSmearLepPt enables smearing of lepton pT.
// Enabled by default.
func WithSmearLepPt(d bool) Option {
	return func(cfg *config) {
		cfg.smearLepPt = d
	}
}

// WithSmearLepTheta enables smearing of lepton
// theta (polar) angle. Enabled by default.
func WithSmearLepTheta(d bool) Option {
	return func(cfg *config) {
		cfg.smearLepTheta = d
	}
}

// WithSmearLepAzimtu enables smearing of lepton
// azimuth angle. Enabled by default.
func WithSmearLepAzimu(d bool) Option {
	return func(cfg *config) {
		cfg.smearLepAzimu = d
	}
}

// WithSmearLepPt enables smearing of jet pT
// (actually jet energy). Enabled by default.
func WithSmearJetPt(d bool) Option {
	return func(cfg *config) {
		cfg.smearJetPt = d
	}
}

// WithSmearJetTheta enables smearing of jet
// theta (polar) angle. Enabled by default.
func WithSmearJetTheta(d bool) Option {
	return func(cfg *config) {
		cfg.smearJetTheta = d
	}
}

// WithSmearJetAzimtu enables smearing of jet
// azimuth angle. Enabled by default.
func WithSmearJetAzimu(d bool) Option {
	return func(cfg *config) {
		cfg.smearJetAzimu = d
	}
}

// WithSmearN sets the number of iteration used in the
// smearing of the kinematics. Default is 10.
func WithSmearN(n int) Option {
	return func(cfg *config) {
		cfg.smearN = n
	}
}

// WithDebug enables debugging print messages.
// Default is false.
func WithDebug(d bool) Option {
	return func(cfg *config) {
		cfg.debug = d
	}
}
