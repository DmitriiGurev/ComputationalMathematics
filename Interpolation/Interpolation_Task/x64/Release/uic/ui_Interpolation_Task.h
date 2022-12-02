/********************************************************************************
** Form generated from reading UI file 'Interpolation_Task.ui'
**
** Created by: Qt User Interface Compiler version 5.12.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_INTERPOLATION_TASK_H
#define UI_INTERPOLATION_TASK_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QSplitter>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QWidget>
#include "qcustomplot.h"

QT_BEGIN_NAMESPACE

class Ui_Interpolation_TaskClass
{
public:
    QAction *actionLoad_File;
    QAction *actionGenerate_random_data;
    QAction *actionInfo;
    QWidget *centralWidget;
    QSplitter *splitter;
    QWidget *widget;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label_2;
    QLabel *label_3;
    QCustomPlot *customPlot;
    QWidget *widget1;
    QHBoxLayout *horizontalLayout;
    QLabel *label;
    QComboBox *comboBox;
    QSpacerItem *horizontalSpacer;
    QMenuBar *menuBar;
    QMenu *menuFile;
    QMenu *menuHelp;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *Interpolation_TaskClass)
    {
        if (Interpolation_TaskClass->objectName().isEmpty())
            Interpolation_TaskClass->setObjectName(QString::fromUtf8("Interpolation_TaskClass"));
        Interpolation_TaskClass->resize(1027, 535);
        actionLoad_File = new QAction(Interpolation_TaskClass);
        actionLoad_File->setObjectName(QString::fromUtf8("actionLoad_File"));
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/Icons/open.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionLoad_File->setIcon(icon);
        actionGenerate_random_data = new QAction(Interpolation_TaskClass);
        actionGenerate_random_data->setObjectName(QString::fromUtf8("actionGenerate_random_data"));
        QIcon icon1;
        icon1.addFile(QString::fromUtf8(":/Icons/die.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionGenerate_random_data->setIcon(icon1);
        actionInfo = new QAction(Interpolation_TaskClass);
        actionInfo->setObjectName(QString::fromUtf8("actionInfo"));
        QIcon icon2;
        icon2.addFile(QString::fromUtf8(":/Icons/info.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionInfo->setIcon(icon2);
        centralWidget = new QWidget(Interpolation_TaskClass);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        splitter = new QSplitter(centralWidget);
        splitter->setObjectName(QString::fromUtf8("splitter"));
        splitter->setGeometry(QRect(10, 10, 1011, 431));
        splitter->setOrientation(Qt::Vertical);
        widget = new QWidget(splitter);
        widget->setObjectName(QString::fromUtf8("widget"));
        horizontalLayout_2 = new QHBoxLayout(widget);
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        horizontalLayout_2->setContentsMargins(0, 0, 0, 0);
        label_2 = new QLabel(widget);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(label_2->sizePolicy().hasHeightForWidth());
        label_2->setSizePolicy(sizePolicy);

        horizontalLayout_2->addWidget(label_2);

        label_3 = new QLabel(widget);
        label_3->setObjectName(QString::fromUtf8("label_3"));

        horizontalLayout_2->addWidget(label_3);

        splitter->addWidget(widget);
        customPlot = new QCustomPlot(splitter);
        customPlot->setObjectName(QString::fromUtf8("customPlot"));
        splitter->addWidget(customPlot);
        widget1 = new QWidget(splitter);
        widget1->setObjectName(QString::fromUtf8("widget1"));
        horizontalLayout = new QHBoxLayout(widget1);
        horizontalLayout->setSpacing(6);
        horizontalLayout->setContentsMargins(11, 11, 11, 11);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        horizontalLayout->setContentsMargins(0, 0, 0, 0);
        label = new QLabel(widget1);
        label->setObjectName(QString::fromUtf8("label"));

        horizontalLayout->addWidget(label);

        comboBox = new QComboBox(widget1);
        comboBox->addItem(QString());
        comboBox->addItem(QString());
        comboBox->addItem(QString());
        comboBox->setObjectName(QString::fromUtf8("comboBox"));

        horizontalLayout->addWidget(comboBox);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer);

        splitter->addWidget(widget1);
        Interpolation_TaskClass->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(Interpolation_TaskClass);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 1027, 21));
        menuFile = new QMenu(menuBar);
        menuFile->setObjectName(QString::fromUtf8("menuFile"));
        menuHelp = new QMenu(menuBar);
        menuHelp->setObjectName(QString::fromUtf8("menuHelp"));
        Interpolation_TaskClass->setMenuBar(menuBar);
        mainToolBar = new QToolBar(Interpolation_TaskClass);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        Interpolation_TaskClass->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(Interpolation_TaskClass);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        Interpolation_TaskClass->setStatusBar(statusBar);

        menuBar->addAction(menuFile->menuAction());
        menuBar->addAction(menuHelp->menuAction());
        menuFile->addAction(actionLoad_File);
        menuFile->addAction(actionGenerate_random_data);
        menuHelp->addAction(actionInfo);
        mainToolBar->addAction(actionLoad_File);
        mainToolBar->addAction(actionGenerate_random_data);
        mainToolBar->addAction(actionInfo);

        retranslateUi(Interpolation_TaskClass);

        QMetaObject::connectSlotsByName(Interpolation_TaskClass);
    } // setupUi

    void retranslateUi(QMainWindow *Interpolation_TaskClass)
    {
        Interpolation_TaskClass->setWindowTitle(QApplication::translate("Interpolation_TaskClass", "Interpolation", nullptr));
        actionLoad_File->setText(QApplication::translate("Interpolation_TaskClass", "Open", nullptr));
        actionGenerate_random_data->setText(QApplication::translate("Interpolation_TaskClass", "Generate random data", nullptr));
        actionInfo->setText(QApplication::translate("Interpolation_TaskClass", "Info", nullptr));
        label_2->setText(QApplication::translate("Interpolation_TaskClass", "Data from:", nullptr));
        label_3->setText(QString());
        label->setText(QApplication::translate("Interpolation_TaskClass", "Interpolation Type", nullptr));
        comboBox->setItemText(0, QApplication::translate("Interpolation_TaskClass", "Linear", nullptr));
        comboBox->setItemText(1, QApplication::translate("Interpolation_TaskClass", "Quadratic", nullptr));
        comboBox->setItemText(2, QApplication::translate("Interpolation_TaskClass", "Cubic", nullptr));

        menuFile->setTitle(QApplication::translate("Interpolation_TaskClass", "File", nullptr));
        menuHelp->setTitle(QApplication::translate("Interpolation_TaskClass", "Help", nullptr));
    } // retranslateUi

};

namespace Ui {
    class Interpolation_TaskClass: public Ui_Interpolation_TaskClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_INTERPOLATION_TASK_H
