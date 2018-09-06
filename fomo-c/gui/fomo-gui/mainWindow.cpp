#include <iostream>
#include <sstream>
#include <chrono>

#include <math.h>

#include <QApplication>
#include <QStyle>
#include <QDesktopWidget>
#include <QGridLayout>
#include <QFormLayout>
#include <QTimer>
#include <QImage>
#include <QLabel>
#include <QGroupBox>
#include <QRadioButton>

#include "mainWindow.h"
#include "ui_mainwindow.h"
#include "FoMo.h"

MainWindow::MainWindow(QApplication &app, QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow) {

	// Initializes this window
	// Caller is still responsible for calling show() when necessary

	// General
	ui->setupUi(this);
	ui->menuBar->hide();
	resize(1280, 720);
	setWindowTitle("FoMo rendering GUI");
	setGeometry(QStyle::alignedRect(Qt::LeftToRight, Qt::AlignCenter, size(), app.desktop()->availableGeometry()));

	// Load config

	// FoMo
	FMO = new FoMo::FoMoObject(3);

	// Top-level layout
	QWidget *centralWidget = new QWidget(this);
	QGridLayout *centralLayout = new QGridLayout(); // No need to assign parent here, setLayout already assigns ownership
	centralWidget->setLayout(centralLayout);
	setCentralWidget(centralWidget);

	// Image
	image = new QImage(512, 512, QImage::Format_ARGB32);
	for(int y = 0; y < 512; y++) {
		for(int x = 0; x < 512; x++) {
			int value = int(255*(time + x/512.0 + y/512.0))%256;
			image->setPixel(x, y, qRgb(value, value, value));
		}
	}
	label = new QLabel(centralWidget);
	label->setPixmap(QPixmap::fromImage(*image));
	centralLayout->addWidget(label, 0, 0, -1, 1, Qt::Alignment(Qt::AlignCenter));

	// Control panel initialization
	QTabWidget *controlPanel = new QTabWidget(centralWidget);
	centralLayout->addWidget(controlPanel, 0, 1);
	QWidget *tab;
	QFormLayout *tabLayout;

	// GoftCube and regular grid
	tab = new QWidget(controlPanel);
	controlPanel->addTab(tab, "Goftcube and regular grid");
	tabLayout = new QFormLayout(); // No need to assign parent here, setLayout already assigns ownership
	tab->setLayout(tabLayout);

	// Log
	log = new Logger(this);
	centralLayout->addWidget(log, 1, 1);

	// Log test
	QTimer *timer = new QTimer(log);
	timer->start(33);
	connect(timer, SIGNAL(timeout()), this, SLOT(PrintTest()));

}

void MainWindow::PrintTest() {
	time += 0.033;
	delete image;
	int width = 257 + int((sin(time*2*3.14159/10)/2 + 0.5)*256)%256;
	int height = 257 + 128 + int((sin(time*2*3.14159/10)/2 + 0.5)*128)%128;
	image = new QImage(width, height, QImage::Format_ARGB32);
	for(int y = 0; y < height; y++) {
		for(int x = 0; x < width; x++) {
			int value = int((time*512 + x + y)/512.0*256)%256;
			image->setPixel(x, y, qRgb(value, value, value));
		}
	}
	auto start = std::chrono::high_resolution_clock::now();
	label->setPixmap(QPixmap::fromImage(*image));
	auto duration = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	log->log(std::to_string(duration));
}

MainWindow::~MainWindow() {
	delete FMO;
	delete image;
	delete ui;
}
