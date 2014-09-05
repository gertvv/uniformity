PACKAGE=uniformity_0.1.tar.gz

$(PACKAGE): uniformity/src/*.c uniformity/R/* uniformity/DESCRIPTION uniformity/NAMESPACE
	rm -f uniformity/src/*.o uniformity/src/*.so
	R CMD build uniformity

.PHONY: install

clean:
	rm -f uniformity/src/*.o uniformity/src/*.so uniformity/src/symbols.rds

install: $(PACKAGE)
	R CMD INSTALL uniformity
