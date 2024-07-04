DIR=$(wildcard ./*/)

all: $(DIR)
	@for d in $^; do make -C $$d ; done

clean: $(DIR)
	@for d in $^; do make -C $$d clean ; done

allclean: $(DIR)
	@for d in $^; do make -C $$d allclean ; done

$(DIRS): FORCE

FORCE: ;