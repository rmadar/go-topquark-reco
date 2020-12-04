package sonn

// Options encodes the various settings to pass to
// the reconstruction.
type Option func(cfg *config)

// config contains all the possible options and their values.
type config struct {
	doSmear bool
	nSmear  int
	debug   bool
}

// newConfig returns a config type with a set of default options.
func newConfig() *config {
	cfg := &config{
		doSmear: true,
		nSmear:  10,
		debug:   false,
	}
	return cfg
}

// WithDebug enables debugging print messages
func WithDebug(d bool) Option {
	return func(cfg *config) {
		cfg.debug = d
	}
}

// WithSmearing enables smearing of kinematics
// Default number of iteration is 10, but can be changed
// using WithNsmearing(n) function.
func WithSmearing(d bool) Option {
	return func(cfg *config) {
		cfg.debug = d
	}
}

// WithNsmearing sets the number of iteration used in the
// smearing of the kinematics.
func WithNsmearing(n int) Option {
	return func(cfg *config) {
		cfg.nSmear = n
	}
}
