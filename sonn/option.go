package sonn

// Options encodes the various settings to pass to
// the reconstruction.
type Option func(cfg *config)

// config contains all the possible options and their values.
type config struct {
	rndseed uint64  
	doSmear bool
	nSmear  int
	debug   bool
}

// newConfig returns a config type with a set of default options.
func newConfig() *config {
	cfg := &config{
		rndseed: 1234,
		doSmear: true,
		nSmear:  10,
		debug:   false,
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

// WithSmearing enables smearing of kinematics
// Default number of iteration is 10, but can be changed
// using WithNsmearing(n) function. Default is true.
func WithSmearing(d bool) Option {
	return func(cfg *config) {
		cfg.debug = d
	}
}

// WithNsmearing sets the number of iteration used in the
// smearing of the kinematics. Default is 10.
func WithNsmearing(n int) Option {
	return func(cfg *config) {
		cfg.nSmear = n
	}
}

// WithDebug enables debugging print messages.
// Default is false.
func WithDebug(d bool) Option {
	return func(cfg *config) {
		cfg.debug = d
	}
}
