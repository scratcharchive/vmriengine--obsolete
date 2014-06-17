TEMPLATE =	app

OBJECTS_DIR = build
DESTDIR = .
TARGET = vmriengine

INCLUDEPATH += .
DEPENDPATH += .
SOURCES += vmriengine.cpp
HEADERS += vmrisimulator.h vmrigpuinterface.h
SOURCES += vmrisimulator.cpp vmrigpuinterface.cpp

INCLUDEPATH += utils
DEPENDPATH += utils
HEADERS += textfile.h qjson.h
SOURCES += textfile.cpp qjson.cpp


LIBS += -lqjson
