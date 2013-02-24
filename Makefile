MAKE			= make

SRCDIR			= src
TESTDIR			= tests

all: laof test
	@echo "DONE"

laof:
	$(MAKE) -C $(SRCDIR)

test:
	$(MAKE) -C $(TESTDIR)

.PHONY: clean
clean:
	$(MAKE) -C $(SRCDIR) clean
	$(MAKE) -C $(TESTDIR) clean
