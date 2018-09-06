#include <QApplication>

#include "mainWindow.h"

int main(int argc, char *argv[]) {

	QApplication app(argc, argv);
	MainWindow window(app);
	window.show();

	return app.exec();

}
