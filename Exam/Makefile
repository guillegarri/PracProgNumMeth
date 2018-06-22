

.PHONY: all
all: regular irregular

.PHONY: regular
regular:
	$(MAKE) -C Regular_grid/

.PHONY: irregular
irregular:
	$(MAKE) -C Irregular_grid/

.PHONY: clean
clean:
	cd Regular_grid/ && $(MAKE) clean
	cd Irregular_grid/ && $(MAKE) clean
