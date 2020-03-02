#!/usr/bin/env python
# -*- coding:utf-8 -*-
from __future__ import division

#-----------------------------------------------------
# File: PaintGeneStructure.py
# Date: 2008-06-30
# Description: 
#-----------------------------------------------------

__version__ = '1.0'
__author__  = 'Wubin Qu <quwubin@gmail.com> @ZCGLAB @BMI @CHINA'
__blog__    = 'http://quwubin.cnblogs.com'
__license__ = 'GPL'

import cairo

def drawUTR(ctx, loci, space, hight, start, length, width, chr_orient):
    ctx.set_source_rgb(0.02, 0.05, 0.67)
    ctx.set_line_width(30)

    left = space + int((int(loci[1]) - start)/length * (width - 2*space))
    right = space + int((int(loci[2]) - start)/length * (width - 2*space))
    if chr_orient == '-':
        left = width - left
        right = width - right
    ctx.move_to(left, hight)
    ctx.line_to(right, hight)
    ctx.stroke()
    return ctx

def drawExon(ctx, loci, space, hight, start, length, width, chr_orient, num, p_num, p_start, p_stop, flag):
    ctx.set_source_rgb(0.51, 0.016, 0.12)
    ctx.set_line_width(30)

    left = space + int((int(loci[1]) - start)/length * (width - 2*space))
    right = space + int((int(loci[2]) - start)/length * (width - 2*space))
    if right - left < 0:
        print 'Error found'
        exit()
    elif right - left == 0:
        right = left + 0.3
    if chr_orient == '-':
        left = width - left
        right = width - right
    ctx.move_to(left, hight)
    ctx.line_to(right, hight)
    ctx.stroke()

    # Exon number
    if flag:
        ctx.set_source_rgb(0, 0, 0)
        ctx.select_font_face("Times", cairo.FONT_SLANT_NORMAL)
        ctx.set_font_size(15)
        tmp_left = left + (right - left) / 2 - 4
        ctx.move_to(tmp_left, hight + 30)
        ctx.show_text(str(num))
        ctx.stroke()

    # For probe
    if str(num) == str(p_num):
        ctx.set_line_width(10)
        ctx.set_source_rgb(0, 0, 0)
        if chr_orient == '+':
            left = left + p_start / length
            right = right - p_stop / length
        else:
            left = left - p_start / length
            right = right + p_stop /length
        hight = hight - 30
        ctx.move_to(left, hight)
        ctx.line_to(right, hight)
        ctx.stroke()

        ctx.set_source_rgb(0, 0, 0)
        ctx.select_font_face("Times", cairo.FONT_SLANT_NORMAL)
        ctx.set_font_size(15)
        text = 'Exon:%s' % num

        length = 0
        for i in range(len(text)):
            length = length + CNW.font20[text[i]]

        text = 'Exon: %s' % num
        tmp_left = left + (right - left) / 2 - length / 2
        ctx.move_to(tmp_left, hight - 10)
        ctx.show_text(text)
        ctx.stroke()


    return ctx

def drawIntro(ctx, a, b, space, hight, start, length, width, chr_orient):
    ctx.set_source_rgb(0, 0, 0)
    ctx.set_line_width(2)
    left = space + int((int(a[2]) - start)/length * (width - 2*space))
    right = space + int((int(b[1]) - start)/length * (width - 2*space))
    if chr_orient == '-':
        left = width - left
        right = width - right
    ctx.move_to(left, hight)
    ctx.line_to(right, hight)
    ctx.stroke()
    return ctx

def forwardTriangle(ctx, pos_x, pos_y):
    ''' Drawing the trangle '''

    ctx.set_source_rgb(0.51, 0.016, 0.12)
    ctx.move_to(pos_x, pos_y)
    ctx.line_to(pos_x - 20, pos_y - 10)
    ctx.rel_line_to(0, 20)
    ctx.close_path()
    ctx.fill()
    ctx.stroke()
    return ctx


def reverseTriangle(ctx, pos_x, pos_y):
    ''' Drawing the trangle '''
    ctx.set_source_rgb(0.51, 0.016, 0.12)
    ctx.move_to(pos_x, pos_y)
    ctx.line_to(pos_x + 20, pos_y - 10)
    ctx.rel_line_to(0, 20)
    ctx.close_path()
    ctx.fill()
    ctx.stroke()
    return ctx

def drawCoordinate(ctx, start, stop, space, width, height):
    ''' Drawing the trangle '''

    # left
    length = 0
    for i in range(len(str(start))):
        length = length + CNW.font20[str(start)[i]]

    pos_x = space - 20 - 10 - length
    pos_y = height + 8

    ctx.move_to(pos_x, pos_y)
    ctx.set_source_rgb(0, 0, 0)
    ctx.select_font_face("Times", cairo.FONT_SLANT_NORMAL)
    ctx.set_font_size(20)
    ctx.show_text(str(start))
    ctx.stroke()

    # right
    pos_x = width - space + 20 + 10
    pos_y = height + 8

    ctx.move_to(pos_x, pos_y)
    ctx.set_source_rgb(0, 0, 0)
    ctx.select_font_face("Times", cairo.FONT_SLANT_NORMAL)
    ctx.set_font_size(20)
    ctx.show_text(str(stop))
    ctx.stroke()

    return ctx


def drawTriangle(ctx, chr_orient, height, space, width, start, stop):
    ''' Drawing the triangle for orientation of chromesome'''
    if chr_orient == '+':
        # left
        pos_x = space
        pos_y = height
        ctx = forwardTriangle(ctx, pos_x, pos_y)
        # right
        pos_x = width - space + 20
        pos_y = height 
        ctx = forwardTriangle(ctx, pos_x, pos_y)
        ctx = drawCoordinate(ctx, start, stop, space, width, height)
    else:
        # left
        pos_x = space - 20
        pos_y = height
        ctx = reverseTriangle(ctx, pos_x, pos_y)
        # right
        pos_x = width - space
        pos_y = height
        ctx = reverseTriangle(ctx, pos_x, pos_y)
        # Exchange the start and stop
        ctx = drawCoordinate(ctx, stop, start, space, width, height)

    return ctx

def drawFootnote(ctx, height, width, space):
    ''' Drawing the foot note below the gene structure pic '''

    ctx.set_source_rgb(0.51, 0.016, 0.12)
    ctx.move_to(space, height)
    ctx.set_line_width(25)
    ctx.line_to(space + 20, height)
    ctx.stroke()

    ctx.move_to(space + 30, height + 12.5)
    ctx.set_source_rgb(0, 0, 0)

    ctx.select_font_face("Times", cairo.FONT_SLANT_NORMAL)
    ctx.set_font_size(20)
    ctx.show_text("- coding region")
    x, y = ctx.get_current_point()
    ctx.stroke()

    ctx.move_to(x + 40, y - 12.5)
    ctx.set_source_rgb(0.02, 0.05, 0.67)
    x, y = ctx.get_current_point()
    ctx.set_line_width(25)
    ctx.line_to(x + 20, y)
    x, y = ctx.get_current_point()
    ctx.stroke()

    ctx.move_to(x + 10, y + 12.5)
    ctx.set_source_rgb(0, 0, 0)

    ctx.show_text("- untranslated region")
    x, y = ctx.get_current_point()
    ctx.stroke()

    ctx.move_to(x + 40, y - 12.5)
    ctx.set_source_rgb(0, 0, 0)
    ctx.set_line_width(2)
    ctx.line_to(x + 60, height)
    ctx.stroke()

    ctx.move_to(x + 70, y)
    ctx.show_text("- Intro")
    ctx.stroke()

    return ctx

def drawTitle(ctx, space, height, width, gene_name):
    length = 0
    for i in range(len(gene_name)):
        length = length + CNW.font30[gene_name[i]]

    start_pos = (width - length) / 2 
    ctx.move_to(start_pos, height)
    ctx.set_source_rgb(0, 0, 0)
    ctx.select_font_face("Times", cairo.FONT_SLANT_NORMAL)
    ctx.set_font_size(30)
    ctx.show_text(gene_name)

    ctx.stroke()


    return ctx

def drawSide(ctx, space, height, width):

    ctx.move_to(space - 30, height + 10)
    ctx.set_source_rgb(0, 0, 0)
    ctx.select_font_face("Times", cairo.FONT_SLANT_NORMAL)
    ctx.set_font_size(30)
    ctx.show_text("5'")

    ctx.move_to(width - space + 10, height + 10)
    ctx.show_text("3'")

    ctx.stroke()


    return ctx

def paint(paint_name, loci_array, start, stop, chr_orient, p_num, p_start, p_stop):
    start = int(start)
    stop = int(stop)

    WIDTH, HEIGHT = 1000, 200
    space = 200

    length = stop - start

    # Setup Cairo
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, HEIGHT)
    ctx = cairo.Context(surface)

    tmp_loci_array = []
    for loci in loci_array:
        if loci[0] == 'UTR':
            print loci
            ctx = drawUTR(ctx, loci, space, HEIGHT/2, start, length, WIDTH, chr_orient)
        else:
            tmp_loci_array.append(loci)

    num = 0
    flag = 0
    for loci in tmp_loci_array:
        num = num + 1
        if num == 1 or num == len(tmp_loci_array):
            flag = 1
        ctx = drawExon(ctx, loci, space, HEIGHT/2, start, length, WIDTH, chr_orient, num, p_num, p_start, p_stop, flag)
        flag = 0

    for i in range(len(loci_array)):
        a = loci_array[i]

        if (i+1) < len(loci_array):
            b = loci_array[i+1]
        else:
            b = ['CDS', stop, stop]

        ctx = drawIntro(ctx, a, b, space, HEIGHT/2, start, length, WIDTH, chr_orient)
        i = i + 1

    ctx = drawTitle(ctx, space, 30, WIDTH, paint_name)
    ctx = drawSide(ctx, space, HEIGHT/2, WIDTH)
    ctx = drawTriangle(ctx, chr_orient, 60, space, WIDTH, start, stop)
    ctx = drawFootnote(ctx, 180, WIDTH, space)
    # Output a PNG file
    paint_name = paint_name + '.png'
    surface.write_to_png(paint_name)


def main():
    paint()

if __name__ == '__main__':
    main()
