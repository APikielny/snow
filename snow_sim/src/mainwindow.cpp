#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "mpm_solver.h"
#include <iostream>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
//    ui->setupUi(this);

    mpm_solver sim;
    float dt = 1e-4;
    sim.initialize();
    while(true) {
        sim.update(dt);
        std::cout << "update" <<std::endl;
    }

}

MainWindow::~MainWindow()
{
    delete ui;
}

