.PHONY: all build

GO := go1.15.6
RM := /bin/rm

build: libcapi.so

clean:
	$(RM) -f ./libcapi.so ./libcapi.h

libcapi.so: capi.go
	$(GO) build -v -o libcapi.so -buildmode=c-shared

test: build test-root.C
	root.exe -b -q ./test-root.C
