#include "Interpolation_Task.h"

#include <QFileDialog>
#include <QString>
#include <QMessageBox>

#include <cstdlib>
#include <iostream>
#include <time.h>
#include <fstream>

#include "interpolationalgorithms.h"

using namespace std;

Interpolation_Task::Interpolation_Task(QWidget *parent)
    : QMainWindow(parent)
{
    ui.setupUi(this);

    setWindowIcon(QIcon(":/Icons/line-graph.png"));

    Interpolation_Task::makePlot();
}

Interpolation_Task::~Interpolation_Task()
{}

void Interpolation_Task::makePlot()
{
    QVector<double> x;
    QVector<double> y;

    if (ui.comboBox->currentText() == "Linear")
    {
        Linear linear;
        linear.LoadData(filename);

        x = QVector<double>::fromStdVector(linear.GetOutput_x());
        y = QVector<double>::fromStdVector(linear.GetOutput_y());
    }

    if (ui.comboBox->currentText() == "Quadratic")
    {
        Quadratic quadratic;
        quadratic.LoadData(filename);

        x = QVector<double>::fromStdVector(quadratic.GetOutput_x());
        y = QVector<double>::fromStdVector(quadratic.GetOutput_y());
    }

    if (ui.comboBox->currentText() == "Cubic")
    {
        Cubic cubic;
        cubic.LoadData(filename);

        x = QVector<double>::fromStdVector(cubic.GetOutput_x());
        y = QVector<double>::fromStdVector(cubic.GetOutput_y());
    }
    ui.customPlot->addGraph();
    ui.customPlot->graph(0)->setData(x, y);
    ui.customPlot->xAxis->setLabel("x");
    ui.customPlot->yAxis->setLabel("y");
    ui.customPlot->xAxis->rescale();
    ui.customPlot->yAxis->rescale();
    ui.customPlot->replot();
}

void Interpolation_Task::on_actionLoad_File_triggered()
{
    filename = QFileDialog::getOpenFileName(
        this,
        tr("Open File"),
        "C://",
        "Text File (*.txt)").toStdString();
    QString f = "<i>" + QFileInfo(QString::fromStdString(filename)).fileName() + "</i>";
    ui.label_3->setText(f);
    makePlot();
}

void Interpolation_Task::on_actionGenerate_random_data_triggered()
{
    srand(time(nullptr));

    ofstream rnd("random.txt");

    if (rnd.is_open())
    {
        double value;
        for (int i = 0; i < 30; i++)
        {
            value = -1 + 2 * rand() / (double)RAND_MAX;

            rnd << i << " " << value << endl;
        }
        rnd.close();
    }
    filename = "random.txt";
    ui.label_3->setText("<i>Random data</i>");

    makePlot();
}

void Interpolation_Task::on_comboBox_currentIndexChanged(const QString& arg1)
{
    makePlot();
}

void Interpolation_Task::on_actionInfo_triggered()
{
    QMessageBox::information(this, tr("Info"),"<h2>Hello</h2> \
        <p> This app offers three polynomial interpolation algorithms: </p> \
        <ul> \
        <li> &nbsp;Linear (&ge;2 points) </li> \
        <li> &nbsp;Quadratic (&ge;3 points) </li> \
        <li> &nbsp;Cubic (&ge;4 points) </li> \
        </ul> \
        <p> The data can be loaded from a .txt file or generated randomly. The generator creates a random set of 30 points with values in [-1, 1].</p> \
        <h2>Important</h2>  \
        <p> Loaded data must meet the following requirements: </p> \
        <ul> \
        <li> &nbsp;The file path contains only ASCII characters</li> \
        <li> &nbsp;<em> h = x[i + 1] - x[i]</em> is the same for all <em>i</em></li> \
        </ul>");
}