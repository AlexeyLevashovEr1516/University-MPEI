;******************************************************************************
;   This file is a basic template for creating relocatable assembly code for  *
;   a PIC18F2520. Copy this file into your project directory and modify or    *
;   add to it as needed.                                                      *
;                                                                             *
;   The PIC18FXXXX architecture allows two interrupt configurations. This     *
;   template code is written for priority interrupt levels and the IPEN bit   *
;   in the RCON register must be set to enable priority levels. If IPEN is    *
;   left in its default zero state, only the interrupt vector at 0x008 will   *
;   be used and the WREG_TEMP, BSR_TEMP and STATUS_TEMP variables will not    *
;   be needed.                                                                *
;                                                                             *
;   Refer to the MPASM User's Guide for additional information on the         *
;   features of the assembler and linker.                                     *
;                                                                             *
;   Refer to the PIC18Fxx20 Data Sheet for additional                         *
;   information on the architecture and instruction set.                      *
;                                                                             *
;******************************************************************************
;                                                                             *
;    Filename:                                                                *
;    Date:                                                                    *
;    File Version:                                                            *
;                                                                             *
;    Author:                                                                  *
;    Company:                                                                 *
;                                                                             * 
;******************************************************************************
;                                                                             *
;    Files required: P18F2520.INC                                             *
;                                                                             *
;******************************************************************************

	LIST P=18F2520, F=INHX32 ;directive to define processor and file format
	#include <p18f2520.INC>	 ;processor specific variable definitions

;******************************************************************************
;Configuration bits
;Microchip has changed the format for defining the configuration bits, please 
;see the .inc file for futher details on notation.  Below are a few examples.



;   Oscillator Selection:
    CONFIG	OSC = INTIO67, WDT = OFF, PBADEN = OFF;, MCLRE = ON, DEBUG = ON             ;LP

;******************************************************************************
;Variable definitions
; These variables are only needed if low priority interrupts are used. 
; More variables may be needed to store other special function registers used
; in the interrupt routines.

		UDATA

WREG_TEMP	RES	1	;variable in RAM for context saving 
STATUS_TEMP	RES	1	;variable in RAM for context saving
BSR_TEMP	RES	1	;variable in RAM for context saving

		UDATA_ACS

LIMIT_GREEN		    RES	    1			    
LIMIT_PULSE_RED		RES		1		    
Voltage				RES		1		

;******************************************************************************
;EEPROM data
; Data to be programmed into the Data EEPROM is defined here

DATA_EEPROM	CODE	0xf00000

		DE	"Test Data",0,1,2,3,4,5

;******************************************************************************
;Reset vector
; This code will start executing when a reset occurs.

RESET_VECTOR	CODE	0x0000

		goto	Main		;go to start of main code

;******************************************************************************
;High priority interrupt vector
; This code will start executing when a high priority interrupt occurs or
; when any interrupt occurs if interrupt priorities are not enabled.

HI_INT_VECTOR	CODE	0x0008

		bra	HighInt		;go to high priority interrupt routine

;******************************************************************************
;Low priority interrupt vector
; This code will start executing when a low priority interrupt occurs.
; This code can be removed if low priority interrupts are not used.

LOW_INT_VECTOR	CODE	0x0018

		bra	LowInt		;go to low priority interrupt routine

;******************************************************************************
;High priority interrupt routine
; The high priority interrupt code is placed here.

		CODE

HighInt:
;	*** high priority interrupt code goes here ***
	RETFIE FAST
	
;******************************************************************************
;Low priority interrupt routine
; The low priority interrupt code is placed here.
; This code can be removed if low priority interrupts are not used.

LowInt:
		movff	STATUS,STATUS_TEMP	;save STATUS register
		movff	WREG,WREG_TEMP		;save working register
		movff	BSR,BSR_TEMP		;save BSR register

;	*** low priority interrupt code goes here ***


		movff	BSR_TEMP,BSR		;restore BSR register
		movff	WREG_TEMP,WREG		;restore working register
		movff	STATUS_TEMP,STATUS	;restore STATUS register
		retfie

;******************************************************************************
;Start of main program
; The main program code is placed here.

Main:
		
	; Конфигурация осциллятора
	MOVLW b'10000000'; F = 8 MHz
	MOVWF OSCTUNE; OSCILLATOR TUNING REGISTER
	MOVLW b'01000011'; Fosc = 1 MHz
	MOVWF OSCCON; OSCILLATOR CONTROL REGISTER
	
	; Конфигурация таймера
	MOVLW b'00000111'; 16-битный, Т=256*4/1MHz=1мс
	MOVWF T0CON; TIMER0 CONTROL REGISTER
	
	; Конфигурация АЦП
	MOVLW b'00010100'; Включаем АЦП
	MOVWF ADCON2
	MOVLW b'00001110'; Аналоговый вход AN0, остальные - цифровые
	MOVWF ADCON1
	MOVLW b'00000001'; Используем канал 0
	MOVWF ADCON0
	
	; Конфигурация прерывания
	BSF INTCON,INT0IE; Нажатие кнопки вызывает прерывания
	BCF INTCON2,INTEDG0; По отрицательному фронту
	BSF RCON,IPEN; Разрешение прерываний
	BSF INTCON,GIEH; Разрешение глобальных прерываний
	
	; Порты входных индикаторов
	BSF TRISA,RA0
	; Порты выходных индикаторов
	BCF TRISB,RB3
	BCF TRISB,RB4
	
	CLRF LATB; Сброс всех значений в регистре LATB
	
	;Указание пределов сравнения для 12В и 14.5В
	MOVLW .153; Запись предела для 12V в WREG
	MOVWF LIMIT_GREEN; Запись из WREG в переменную LIMIT_GREEN
	MOVLW .185; Запись передела для 14.5V в WREG
	MOVWF LIMIT_PULSE_RED; Запись из WREG в переменную LIMIT_PULSE_RED

M1:	
		BSF ADCON0,1; Запуск АЦП
MAD:	BTFSS ADCON0,1; Проверка на завершение оцифровки
		    BRA MAD
		MOVFF ADRESH,Voltage; Запись полученного значения в Voltage 
	
	MOVF LIMIT_GREEN,W; Запись LIMIT_GREEN в WREG
	CPFSLT Voltage; Сравнение Voltage и WREG
	    BRA M2; Идти на метку М2, если Voltage > WREG
	BCF LATB,RB3; Выключить RB3
	BSF LATB,RB4; Включить RB4
	BRA M1; Вернуться на метку М1

M2: MOVF LIMIT_PULSE_RED,W; Запись LIMIT_PULSE_RED в WREG
	CPFSLT Voltage; Сравнение Voltage и WREG
	    BRA M3; Идти на метку М2, если Voltage > WREG
	BSF LATB,RB3; Включить RB3
	BCF LATB,RB4; Выключить RB4
	BRA M1; Вернуться на метку М1

	; Так как период таймера Т=1мс, а для переключения светодиода необходимо 500 мс, 
	; значит количество тактов таймера будет 500, следовательно начальное значение таймера 
	; 65535 - 500 = 65035, что будет 0xFE0B
	
M3: BCF LATB,RB3; Выключить RB3
	MOVLW 0xFE; Запись 0xF0 в WREG
	MOVWF TMR1H; Запись из WREG в перемнную TMR1H
	MOVLW 0x0B; Запись 0x0B в WREG
	MOVWF TMR1L; Запись из WREG в перемнную TMR1L
	BCF PIR1,TMR1IF; Сброс флага таймера
	BSF INT0CON,TMR0ON; Запустить таймер
M4:	BTFSS PIR1,TMR0IF; Проверка на переполнение таймера (500 мс)
		BRA M4; Вернуться на метку М4
	BCF INTCON,TMR0IF; Сброс флага переполнения таймера
	BTG LATB,RB4; Инверсировать бит RB4
	BRA M1; Вернуться на метку М1

;******************************************************************************
;End of program
		
    END


