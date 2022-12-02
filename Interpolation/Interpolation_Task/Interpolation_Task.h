#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_Interpolation_Task.h"

#include <string>

class Interpolation_Task : public QMainWindow
{
    Q_OBJECT

public:
    Interpolation_Task(QWidget *parent = nullptr);
    ~Interpolation_Task();

protected:
    std::string filename = "";

private slots:
    void makePlot();

    void on_actionLoad_File_triggered();

    void on_actionGenerate_random_data_triggered();

    void on_comboBox_currentIndexChanged(const QString& arg1);

    void on_actionInfo_triggered();

private:
    Ui::Interpolation_TaskClass ui;
};
